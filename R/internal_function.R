# 内部函数
.cal_prob <- function(expr_matrix, gs, wt_gs, 
                      pan_scores) {
  n_cells <- ncol(expr_matrix)
  n_types <- length(gs)
  cell_ids <- colnames(expr_matrix)
  
  # 初始化概率矩阵
  prob_matrix <- matrix(1/n_types, nrow = n_cells, ncol = n_types,
                        dimnames = list(cell_ids, names(gs)))
  
  # 处理marker权重
  calculate_scientific_weights <- function(weight_levels) {
    weight_values <- sapply(weight_levels, function(level) {
      if (is.numeric(level)) return(level)
      switch(tolower(level),
             "high" = 1.0,
             "medium" = 0.7,
             "low" = 0.3,
             1.0)
    })
    weight_values / mean(weight_values)
  }
  
  # 创建权重矩阵 (仅使用gs中的基因)
  all_marker_genes <- unique(unlist(gs))
  weight_matrix <- matrix(0, nrow = n_types, ncol = length(all_marker_genes),
                          dimnames = list(names(gs), all_marker_genes))
  
  # 步骤1: 初始化所有基础标记基因的权重为low (0.3)
  for (type in names(gs)) {
    weight_matrix[type, gs[[type]]] <- 0.3
  }
  
  # 步骤2: 应用用户提供的权重 (覆盖默认值)
  if (!is.null(wt_gs) && length(wt_gs) > 0) {
    for (type in names(wt_gs)) {
      if (type %in% names(gs)) {
        weights <- wt_gs[[type]]
        sci_weights <- calculate_scientific_weights(weights)
        
        # 仅处理存在于基础标记基因中的基因
        valid_genes <- intersect(names(weights), gs[[type]])
        if (length(valid_genes) > 0) {
          weight_matrix[type, valid_genes] <- sci_weights[valid_genes]
        }
      }
    }
  }
  
  # 计算细胞类型概率
  for (i in seq_len(n_types)) {
    type <- names(gs)[i]
    relevant_genes <- gs[[type]]
    valid_genes <- intersect(relevant_genes, rownames(expr_matrix))
    
    if (length(valid_genes) > 0) {
      weights <- weight_matrix[type, valid_genes]
      expr_vals <- expr_matrix[valid_genes, , drop = FALSE]
      expr_vals[expr_vals < 0] <- 0
      
      weighted_geo_mean <- exp(base::colSums(weights * log1p(expr_vals)) / sum(weights))
      prob_matrix[, i] <- weighted_geo_mean
    } else {
      prob_matrix[, i] <- 0
    }
  }
  
  # Softmax归一化
  prob_matrix <- exp(prob_matrix - matrixStats::rowMaxs(prob_matrix))
  prob_matrix <- prob_matrix / base::rowSums(prob_matrix)
  
  # 生成结果
  results <- data.frame(
    cell_id = cell_ids,
    predicted.label = colnames(prob_matrix)[max.col(prob_matrix)],
    predicted.score = matrixStats::rowMaxs(prob_matrix),
    stringsAsFactors = FALSE
  )
  
  # 添加概率矩阵
  results <- cbind(results, as.data.frame(prob_matrix))
  
  return(results)
}

.extract_gs <- function(gs, wt_gs) {
  marker_info <- list()
  
  for (celltype in names(gs)) {
    # 获取基础标记基因
    base_genes <- gs[[celltype]]
    
    # 如果有权重信息，获取加权基因
    weighted_genes <- if (!is.null(wt_gs) && celltype %in% names(wt_gs)) {
      names(wt_gs[[celltype]])
    } else {
      NULL
    }
    
    # 确定最终使用的标记基因
    if (!is.null(weighted_genes)) {
      # 取基础基因和加权基因的交集
      final_genes <- intersect(base_genes, weighted_genes)
      if (length(final_genes) == 0) {
        warning("No overlapping marker genes for celltype: ", celltype,
                ". Using all base markers.")
        final_genes <- base_genes
      }
    } else {
      final_genes <- base_genes
    }
    
    marker_info[[celltype]] <- final_genes
  }
  
  return(marker_info)
}