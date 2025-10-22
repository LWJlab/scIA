#' Single-cell re-identification analysis (single)
#' This function performs single-cell re-identification by evaluating 
#' whether cells expressing a sufficient number of marker genes should 
#' be reclassified to a target cell type.
#' @param object Seurat object
#' @param cell_ids Name of cells to be analyzed
#' @param current_labels The current cell type of the cells
#' @param compared_labels The cell type to be identified
#' @param gs List of marker genes
#' @param threshold The number of marker genes
#' @param add_gs_exp Whether to print the gene expression in the results
#' @return A dataframe containing the final results
#' @examples
#' \dontrun{
#' scReType_single(sce, 
#' cell_ids = cell_ids, 
#' current_labels = 'Undetermined',
#' compared_labels = 'Neutrophils',
#' gs = gs, 
#' threshold = 2,
#' add_gs_exp = FALSE
#' }
#' @importFrom Seurat GetAssayData
#' @importFrom Matrix colSums
#' @export
scReType_single <- function(object, 
                            cell_ids, 
                            current_labels, 
                            compared_labels, 
                            gs, 
                            threshold = 3,
                            add_gs_exp = FALSE) {
  
  # 添加细胞数量提示信息
  message("Processing ", length(cell_ids), " cells for multi cell type re-identification")
  
  # 验证输入参数
  if (length(compared_labels) != 1) {
    stop("compared_labels must be a single cell type for scReType_single")
  }
  
  if (length(current_labels) == 1) {
    current_labels <- rep(current_labels, length(cell_ids))
    message("Expanded single current_label to all cells: '", current_labels[1], "'")
  }
  
  if (length(cell_ids) != length(current_labels)) {
    stop("cell_ids (length ", length(cell_ids), 
         ") and current_labels (length ", length(current_labels), 
         ") must have the same length")
  }
  
  if (!all(cell_ids %in% colnames(object))) {
    missing_cells <- setdiff(cell_ids, colnames(object))
    stop("Some cell_ids are not found in the Seurat object: ", 
         paste(head(missing_cells, 5), collapse = ", "),
         ifelse(length(missing_cells) > 5, "...", ""))
  }
  
  # 获取表达矩阵
  expr_matrix <- Seurat::GetAssayData(object, slot = "data")[, cell_ids, drop = FALSE]
  
  # 初始化结果数据框
  result_df <- data.frame(
    Input_confidence = threshold,
    predicted.label = current_labels,
    stringsAsFactors = FALSE
  )
  
  # 添加置信度列
  conf_colname <- paste0("Confidence_", compared_labels)
  result_df[[conf_colname]] <- 0
  
  # 确保标记基因在表达矩阵中存在
  valid_markers <- intersect(gs, rownames(expr_matrix))
  if (length(valid_markers) == 0) {
    warning("None of the provided markers are found in the expression matrix")
    return(result_df)
  }
  
  message("Processing '", compared_labels, "' with ", length(valid_markers), 
          " valid markers: ", paste(valid_markers, collapse = ", "))
  
  # 计算每个细胞表达标记基因的数量
  marker_expr <- expr_matrix[valid_markers, , drop = FALSE] > 0
  expr_count <- Matrix::colSums(marker_expr)
  
  # 存储置信度
  result_df[[conf_colname]] <- expr_count
  
  # 重新分类细胞
  reclassify_idx <- expr_count >= threshold
  if (any(reclassify_idx)) {
    result_df$predicted.label[reclassify_idx] <- compared_labels
    message("Reclassified ", sum(reclassify_idx), " cells to '", compared_labels, "'")
  }
  
  # 如果用户选择添加标记基因表达量
  if (add_gs_exp && length(valid_markers) > 0) {
    # 获取这些标记基因的表达量
    marker_expr_data <- as.data.frame(t(as.matrix(expr_matrix[valid_markers, , drop = FALSE])))
    
    # 为每个标记基因添加细胞类型前缀
    colnames(marker_expr_data) <- paste0(compared_labels, "_", colnames(marker_expr_data))
    
    # 将表达量添加到结果数据框
    result_df <- cbind(result_df, marker_expr_data)
    
    message("Added expression data for ", length(valid_markers), 
            " markers of '", compared_labels, "' type")
  }
  
  # 确保行顺序与输入cell_ids一致
  rownames(result_df) <- cell_ids
  
  return(result_df)
}







#' Single-cell re-identification analysis (multi)
#' This function performs multi-cell re-identification by evaluating 
#' cells against multiple candidate cell types simultaneously and assigning 
#' the most probable label based on marker gene expression patterns.
#' @param object Seurat object
#' @param cell_ids Name of cells to be analyzed
#' @param current_labels The current cell type of the cells
#' @param compared_labels The cell type to be identified
#' @param gs List of marker genes
#' @param shared_gs Shared marker genes for cell types to be identified
#' @param add_gs_exp Whether to print the gene expression in the results
#' @return A dataframe containing the final results
#' @examples
#' \dontrun{
#' scReType_multi(sce, 
#' cell_ids = cell_ids, 
#' current_labels = 'Undetermined',
#' compared_labels = 'Neutrophils',
#' gs = gs, 
#' shared_gs = shared_gs,
#' add_gs_exp = FALSE
#' }
#' @importFrom Seurat GetAssayData
#' @importFrom Matrix colSums
#' @export
scReType_multi <- function(object, 
                           cell_ids, 
                           current_labels, 
                           compared_labels, 
                           gs, 
                           shared_gs = NULL,
                           add_gs_exp = FALSE) {
  # 添加细胞数量提示信息
  message("Processing ", length(cell_ids), " cells for multi cell type re-identification")
  
  # 验证输入参数
  if (length(compared_labels) < 2) {
    stop("scReType_multi requires at least 2 cell types for comparison")
  }
  
  if (length(compared_labels) != length(gs)) {
    stop("compared_labels and gs must have the same length")
  }
  
  if (!all(names(gs) %in% compared_labels)) {
    stop("All names in gs must be present in compared_labels")
  }
  
  # 验证shared_gs
  if (!is.null(shared_gs)) {
    if (!is.list(shared_gs)) {
      stop("shared_gs must be a list")
    }
    
    # 检查所有共享标记基因对应的细胞类型都在compared_labels中
    all_shared_celltypes <- unique(unlist(shared_gs))
    if (!all(all_shared_celltypes %in% compared_labels)) {
      missing_types <- setdiff(all_shared_celltypes, compared_labels)
      stop("Some cell types in shared_gs are not in compared_labels: ", 
           paste(missing_types, collapse = ", "))
    }
  }
  
  if (length(current_labels) == 1) {
    current_labels <- rep(current_labels, length(cell_ids))
    message("Expanded single current_label to all cells: '", current_labels[1], "'")
  }
  
  if (length(cell_ids) != length(current_labels)) {
    stop("cell_ids (length ", length(cell_ids), 
         ") and current_labels (length ", length(current_labels), 
         ") must have the same length")
  }
  
  if (!all(cell_ids %in% colnames(object))) {
    missing_cells <- setdiff(cell_ids, colnames(object))
    stop("Some cell_ids are not found in the Seurat object: ", 
         paste(head(missing_cells, 5), collapse = ", "),
         ifelse(length(missing_cells) > 5, "...", ""))
  }
  
  # 获取表达矩阵
  expr_matrix <- Seurat::GetAssayData(object, slot = "data")[, cell_ids, drop = FALSE]
  
  # 初始化结果数据框
  result_df <- data.frame(
    predicted.label = current_labels,
    stringsAsFactors = FALSE
  )
  
  # 为每种细胞类型添加置信度列
  for (celltype in compared_labels) {
    conf_colname <- paste0("Confidence_", celltype)
    result_df[[conf_colname]] <- 0
  }
  
  # 创建一个矩阵存储每种细胞类型的置信度
  confidence_matrix <- matrix(0, 
                              nrow = length(cell_ids), 
                              ncol = length(compared_labels),
                              dimnames = list(cell_ids, compared_labels))
  
  # 处理共享标记基因 - 为每种细胞类型创建完整的标记基因列表
  full_markers_list <- gs
  
  if (!is.null(shared_gs)) {
    # 为每种细胞类型添加共享标记基因
    for (celltype in compared_labels) {
      # 找到包含该细胞类型的所有共享标记基因
      shared_genes <- names(shared_gs)[sapply(shared_gs, function(x) celltype %in% x)]
      
      # 将这些共享标记基因添加到该细胞类型的标记列表中
      full_markers_list[[celltype]] <- c(full_markers_list[[celltype]], shared_genes)
    }
  }
  
  # 用于存储所有要添加的表达式数据（避免重复添加共享标记基因）
  all_expr_data <- list()
  
  # 对每种细胞类型进行处理
  for (celltype in compared_labels) {
    markers <- full_markers_list[[celltype]]
    
    # 分离普通标记基因和共享标记基因
    regular_markers <- markers[!markers %in% names(shared_gs)]
    shared_markers <- markers[markers %in% names(shared_gs)]
    
    # 处理普通标记基因
    valid_regular_markers <- intersect(regular_markers, rownames(expr_matrix))
    
    # 处理共享标记基因
    valid_shared_markers <- intersect(shared_markers, rownames(expr_matrix))
    
    # 合并所有有效标记基因
    valid_markers <- c(valid_regular_markers, valid_shared_markers)
    
    if (length(valid_markers) == 0) {
      warning("None of the provided markers for '", celltype, "' are found in the expression matrix")
      next
    }
    
    message("Processing '", celltype, "' with ", length(valid_markers), 
            " valid markers: ", paste(valid_markers, collapse = ", "))
    
    # 计算每个细胞表达该类型标记基因的数量
    marker_expr <- expr_matrix[valid_markers, , drop = FALSE] > 0
    expr_count <- Matrix::colSums(marker_expr)
    
    # 存储该细胞类型的置信度
    conf_colname <- paste0("Confidence_", celltype)
    result_df[[conf_colname]] <- expr_count
    confidence_matrix[, celltype] <- expr_count
    
    # 如果用户选择添加标记基因表达量
    if (add_gs_exp) {
      # 处理普通标记基因
      if (length(valid_regular_markers) > 0) {
        # 获取这些标记基因的表达量
        marker_expr_data <- as.data.frame(t(as.matrix(expr_matrix[valid_regular_markers, , drop = FALSE])))
        
        # 为每个标记基因添加细胞类型前缀
        colnames(marker_expr_data) <- paste0(celltype, "_", colnames(marker_expr_data))
        
        # 将表达量添加到所有表达式数据列表
        all_expr_data <- c(all_expr_data, list(marker_expr_data))
        
        message("Added expression data for ", length(valid_regular_markers), 
                " regular markers of '", celltype, "' type")
      }
      
      # 处理共享标记基因（只在第一次遇到时添加）
      if (length(valid_shared_markers) > 0) {
        # 只添加尚未处理的共享标记基因
        new_shared_markers <- setdiff(valid_shared_markers, names(all_expr_data))
        
        if (length(new_shared_markers) > 0) {
          # 获取这些共享标记基因的表达量
          shared_expr_data <- as.data.frame(t(as.matrix(expr_matrix[new_shared_markers, , drop = FALSE])))
          
          # 为每个共享标记基因添加适当的前缀
          for (marker in colnames(shared_expr_data)) {
            # 共享标记基因使用共享细胞类型前缀
            shared_celltypes <- sort(shared_gs[[marker]])
            prefix <- paste0(shared_celltypes, collapse = "_")
            new_colname <- paste0(prefix, "_", marker)
            
            # 重命名列
            colnames(shared_expr_data)[colnames(shared_expr_data) == marker] <- new_colname
            
            # 标记这个共享基因已经处理过
            all_expr_data[[marker]] <- TRUE
          }
          
          # 将共享标记基因表达量添加到所有表达式数据列表
          all_expr_data <- c(all_expr_data, list(shared_expr_data))
          
          message("Added expression data for ", length(new_shared_markers), 
                  " shared markers: ", paste(new_shared_markers, collapse = ", "))
        }
      }
    }
  }
  
  # 将所有表达式数据添加到结果数据框
  if (add_gs_exp && length(all_expr_data) > 0) {
    # 过滤掉非数据框元素（标记）
    expr_data_frames <- all_expr_data[sapply(all_expr_data, is.data.frame)]
    
    if (length(expr_data_frames) > 0) {
      # 合并所有表达式数据
      combined_expr_data <- do.call(cbind, expr_data_frames)
      
      # 将表达式数据添加到结果数据框
      result_df <- cbind(result_df, combined_expr_data)
      
      message("Added expression data for ", ncol(combined_expr_data), " markers in total")
    }
  }
  
  # 确定每个细胞的最高置信度类型
  for (i in seq_along(cell_ids)) {
    conf_scores <- confidence_matrix[i, ]
    max_conf <- max(conf_scores)
    
    if (max_conf > 0) {
      # 找到所有达到最大置信度的类型
      max_types <- names(conf_scores)[conf_scores == max_conf]
      
      if (length(max_types) == 1) {
        # 如果只有一个类型达到最大置信度，则重新分类
        result_df$predicted.label[i] <- max_types
      } else {
        # 如果有多个类型达到相同置信度，保留原标签
        # 可以在这里添加更复杂的决策逻辑
      }
    }
  }
  
  # 确保行顺序与输入cell_ids一致
  rownames(result_df) <- cell_ids
  
  # 统计重新分类的细胞数量
  reclassified <- sum(result_df$predicted.label != current_labels)
  if (reclassified > 0) {
    message("Reclassified ", reclassified, " cells (", 
            round(reclassified/length(cell_ids)*100, 1), "%)")
  }
  
  return(result_df)
}
