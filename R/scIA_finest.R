#' Single-cell intelligent annotation (finest)
#' 
#' This function determines more accurate cell type through a complex decision system 
#' based on the semi-supervised and unsupervised results of scIA function.
#' @param df A dataframe containing the semi-supervised and unsupervised results of scIA function
#' @param pan_mapping Mapping cells to main cell types
#' @return A dataframe containing the final results
#' @examples
#' \dontrun{
#' scIA_finest(df,
#' pan_mapping = pan_mapping
#' }
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
scIA_finest <- function(df,
                        pan_mapping = NULL) {
  # 显示开始信息
  message("Starting scIA_finest: Automated cell type annotation based on scIA default results")
  
  # 定义默认的泛大类细胞类型映射
  default_pan_mapping <- list(
    Pan_epithelial = c("Basal resting", "Suprabasal", "Hillock-like", "Deuterosomal", 
                       "Multiciliated (non-nasal)", "Club (non-nasal)", "Goblet (bronchial)", 
                       "Goblet (subsegmental)", "SMG serous (bronchial)", "SMG mucous",
                       "SMG duct", "Tuft", "Ionocyte", "Neuroendocrine", "AT1", "AT2"),
    
    Pan_endothelial = c("EC general capillary", "EC aerocyte capillary", "EC arterial", 
                        "EC venous pulmonary", "EC venous systemic", "Lymphatic EC mature",
                        "Lymphatic EC differentiating"),
    
    Pan_immune = c("Alveolar macrophages", "Monocyte-derived Mph", "Interstitial Mph perivascular",
                   "Classical monocytes", "Non-classical monocytes", "DC1", "DC2", "Migratory DCs",
                   "Plasmacytoid DCs", "Mast cells", "B cells", "Plasma cells", "CD4 T cells",
                   "CD8 T cells", "NK cells", "Hematopoietic stem cells"),
    
    Pan_stroma = c("Alveolar fibroblasts", "Adventitial fibroblasts", "Peribronchial fibroblasts",
                   "Subpleural fibroblasts", "Myofibroblasts", "Pericytes", "Smooth muscle", "Mesothelium")
  )
  
  # 使用用户提供的映射或默认映射
  if (is.null(pan_mapping)) {
    pan_mapping <- default_pan_mapping
    message("Using default pan-category mapping")
  } else {
    # 检查用户提供的映射是否包含必要的泛大类
    required_categories <- c("Pan_epithelial", "Pan_endothelial", "Pan_immune", "Pan_stroma")
    missing_categories <- setdiff(required_categories, names(pan_mapping))
    
    if (length(missing_categories) > 0) {
      warning("User-provided pan_mapping is missing the following categories: ", 
              paste(missing_categories, collapse = ", "), 
              ". Using default mapping for missing categories.")
      
      # 使用默认值补充缺失的类别
      for (category in missing_categories) {
        if (!category %in% names(pan_mapping)) {
          pan_mapping[[category]] <- default_pan_mapping[[category]]
        }
      }
    }
    message("Using user-provided pan-category mapping")
  }
  
  # 提取泛大类细胞类型列表
  Pan_epithelial <- pan_mapping$Pan_epithelial
  Pan_endothelial <- pan_mapping$Pan_endothelial
  Pan_immune <- pan_mapping$Pan_immune
  Pan_stroma <- pan_mapping$Pan_stroma
  
  # 辅助函数：获取细胞类型所属的泛大类
  get_pan_category <- function(cell_type) {
    if (is.na(cell_type)) return(NA)
    if (cell_type %in% Pan_epithelial) return("Pan_epithelial")
    if (cell_type %in% Pan_endothelial) return("Pan_endothelial")
    if (cell_type %in% Pan_immune) return("Pan_immune")
    if (cell_type %in% Pan_stroma) return("Pan_stroma")
    return(NA)
  }
  
  # 辅助函数：判断是否为特殊情况（AT2与基质，或内皮与基质）
  is_special_case <- function(type1, type2) {
    if (is.na(type1) || is.na(type2)) return(FALSE)
    
    pan1 <- get_pan_category(type1)
    pan2 <- get_pan_category(type2)
    
    if (is.na(pan1) || is.na(pan2)) return(FALSE)
    
    # 上皮与任意基质类型
    if (pan1 == "Pan_epithelial" && pan2 == "Pan_stroma") return(TRUE)
    if (pan2 == "Pan_epithelial" && pan1 == "Pan_stroma") return(TRUE)
    
    # 内皮与任意基质类型
    if (pan1 == "Pan_endothelial" && pan2 == "Pan_stroma") return(TRUE)
    if (pan2 == "Pan_endothelial" && pan1 == "Pan_stroma") return(TRUE)
    
    return(FALSE)
  }
  
  # 辅助函数：处理半监督预测（规则一）
  process_semi_supervised <- function(row) {
    semi_algorithms <- c("COSG_semi", "FindAllMarkers_semi", "GenesorteR_semi", "presto_semi")
    semi_confidences <- c("Conf_COSG_semi", "Conf_FindAllMarkers_semi", "Conf_GenesorteR_semi", "Conf_presto_semi")
    
    predictions <- unlist(row[semi_algorithms])
    confidences <- as.numeric(unlist(row[semi_confidences]))
    
    # 处理NA值
    predictions[is.na(predictions)] <- "Unknown"
    confidences[is.na(confidences)] <- 0
    
    # 统计每种预测出现的次数
    pred_table <- table(predictions)
    max_count <- max(pred_table)
    
    # 规则1.1：三种及以上算法预测相同
    if (max_count >= 3) {
      return(names(pred_table)[which.max(pred_table)])
    }
    
    # 规则1.2：两种算法预测相同，另外两种预测另一种相同类型
    if (max_count == 2 && length(pred_table) == 2) {
      pred1 <- names(pred_table)[1]
      pred2 <- names(pred_table)[2]
      
      # 获取两种预测的置信度（取两种算法中较高的置信度）
      conf1 <- max(confidences[predictions == pred1], na.rm = TRUE)
      conf2 <- max(confidences[predictions == pred2], na.rm = TRUE)
      
      if (conf1 > conf2) return(pred1)
      if (conf2 > conf1) return(pred2)
      return("Undetermined")
    }
    
    # 规则1.3：两种算法预测相同，另外两种预测不同类型
    if (max_count == 2 && length(pred_table) == 3) {
      main_pred <- names(which.max(pred_table))
      other_preds <- names(pred_table)[names(pred_table) != main_pred]
      
      main_conf <- max(confidences[predictions == main_pred], na.rm = TRUE)
      other_confs <- sapply(other_preds, function(p) {
        max(confidences[predictions == p], na.rm = TRUE)
      })
      
      # 如果其他两种预测的置信度都大于0或任意一种等于或大于主要预测的置信度
      if (all(other_confs > 0) || any(other_confs >= main_conf)) {
        return("Low quality")
      } else {
        return(main_pred)
      }
    }
    
    # 规则1.4：四种算法预测四种不同类型
    if (length(unique(predictions)) == 4) {
      return("Low quality")
    }
    
    return("Undetermined") # 默认情况
  }
  
  # 辅助函数：处理无监督预测（规则二）
  process_unsupervised <- function(row) {
    all_algorithms <- c("COSG_All", "FindAllMarkers_All", "GenesorteR_All", "presto_All")
    all_confidences <- c("Conf_COSG_All", "Conf_FindAllMarkers_All", "Conf_GenesorteR_All", "Conf_presto_All")
    
    predictions <- unlist(row[all_algorithms])
    confidences <- as.numeric(unlist(row[all_confidences]))
    
    # 处理NA值
    predictions[is.na(predictions)] <- "Unknown"
    confidences[is.na(confidences)] <- 0
    
    # 统计每种预测出现的次数
    pred_table <- table(predictions)
    max_count <- max(pred_table)
    
    # 规则2.1：三种及以上算法预测相同
    if (max_count >= 3) {
      return(names(pred_table)[which.max(pred_table)])
    }
    
    # 规则2.2：两种算法预测相同，另外两种预测另一种类型
    if (max_count == 2 && length(pred_table) == 2) {
      pred1 <- names(pred_table)[1]
      pred2 <- names(pred_table)[2]
      pan1 <- get_pan_category(pred1)
      pan2 <- get_pan_category(pred2)
      
      conf1 <- max(confidences[predictions == pred1], na.rm = TRUE)
      conf2 <- max(confidences[predictions == pred2], na.rm = TRUE)
      
      # 相同大类
      if (!is.na(pan1) && !is.na(pan2) && pan1 == pan2) {
        if (conf1 > conf2) return(pred1)
        if (conf2 > conf1) return(pred2)
        return("Undetermined")
      }
      
      # 不同大类
      if (conf1 > 0 && conf2 == 0) return(pred1)
      if (conf2 > 0 && conf1 == 0) return(pred2)
      if (conf1 > 0 && conf2 > 0) {
        if (is_special_case(pred1, pred2)) {
          if (conf1 > conf2) return(pred1)
          if (conf2 > conf1) return(pred2)
          return("Undetermined")
        }
        return("Low quality")
      }
      return("Low quality")
    }
    
    # 规则2.3：两种算法预测相同，另外两种预测不同类型
    if (max_count == 2 && length(pred_table) == 3) {
      main_pred <- names(which.max(pred_table))
      other_preds <- names(pred_table)[names(pred_table) != main_pred]
      
      main_conf <- max(confidences[predictions == main_pred], na.rm = TRUE)
      other_confs <- sapply(other_preds, function(p) {
        max(confidences[predictions == p], na.rm = TRUE)
      })
      
      # 如果其他两种预测的置信度都大于0或任意一种等于或大于主要预测的置信度
      if (all(other_confs > 0) || any(other_confs >= main_conf)) {
        return("Low quality")
      } else {
        return(main_pred)
      }
    }
    
    # 规则2.4：四种算法预测四种不同类型
    if (length(unique(predictions)) == 4) {
      return("Low quality")
    }
    
    return("Undetermined") # 默认情况
  }
  
  # 辅助函数：获取最佳细胞类型（四种算法中出现次数最多的类型，至少出现两次且置信度最大）
  get_best_prediction <- function(predictions, confidences) {
    # 处理NA值
    predictions[is.na(predictions)] <- "Unknown"
    confidences[is.na(confidences)] <- 0
    
    pred_table <- table(predictions)
    max_count <- max(pred_table)
    
    # 如果最大出现次数小于2，则没有最佳预测
    if (max_count < 2) {
      return(NA)
    }
    
    # 如果有多个类型出现次数相同，选择置信度最高的
    if (sum(pred_table == max_count) > 1) {
      candidates <- names(pred_table)[pred_table == max_count]
      candidate_confs <- sapply(candidates, function(candidate) {
        max(confidences[predictions == candidate], na.rm = TRUE)
      })
      return(candidates[which.max(candidate_confs)])
    }
    
    return(names(which.max(pred_table)))
  }
  
  # 辅助函数：获取最佳细胞类型的置信度
  get_best_confidence <- function(predictions, confidences, best_pred) {
    # 处理NA值
    predictions[is.na(predictions)] <- "Unknown"
    confidences[is.na(confidences)] <- 0
    
    # 获取最佳预测的置信度（取该类型中所有算法置信度的最大值）
    return(max(confidences[predictions == best_pred], na.rm = TRUE))
  }
  
  # 创建进度条
  message("Processing cells: ")
  pb <- utils::txtProgressBar(min = 0, max = nrow(df), style = 3)
  
  # 主循环：处理每一行数据
  predicted_labels <- vector("character", nrow(df))
  
  for (i in 1:nrow(df)) {
    # 更新进度条
    utils::setTxtProgressBar(pb, i)
    
    row <- df[i, ]
    
    # 获取泛评分并处理NA值
    pan_scores <- c(
      Pan_epithelial = as.numeric(row["Pan_epithelial"]),
      Pan_endothelial = as.numeric(row["Pan_endothelial"]),
      Pan_immune = as.numeric(row["Pan_immune"]),
      Pan_stroma = as.numeric(row["Pan_stroma"])
    )
    pan_scores[is.na(pan_scores)] <- 0
    
    # 获取半监督预测
    semi_algorithms <- c("COSG_semi", "FindAllMarkers_semi", "GenesorteR_semi", "presto_semi")
    semi_predictions <- unlist(row[semi_algorithms])
    
    # 处理NA值
    semi_predictions[is.na(semi_predictions)] <- "Unknown"
    
    # 确定半监督预测的泛大类
    semi_pan_categories <- sapply(semi_predictions, get_pan_category)
    
    # 检查是否所有半监督预测属于同一泛大类
    if (length(unique(na.omit(semi_pan_categories))) == 1 && !all(is.na(semi_pan_categories))) {
      pan_category <- unique(na.omit(semi_pan_categories))[1]
      
      # 规则一：匹配情况
      if (pan_scores[pan_category] > 0 && all(pan_scores[names(pan_scores) != pan_category] == 0)) {
        predicted_labels[i] <- process_semi_supervised(row)
        next
      }
      
      # 规则二：不匹配情况
      if (pan_scores[pan_category] == 0) {
        predicted_labels[i] <- process_unsupervised(row)
        next
      }
    }
    
    # 规则三：不完全匹配或所有评分为0
    # 获取半监督和无监督的最佳预测
    semi_algorithms <- c("COSG_semi", "FindAllMarkers_semi", "GenesorteR_semi", "presto_semi")
    semi_confidences <- c("Conf_COSG_semi", "Conf_FindAllMarkers_semi", "Conf_GenesorteR_semi", "Conf_presto_semi")
    semi_preds <- unlist(row[semi_algorithms])
    semi_confs <- as.numeric(unlist(row[semi_confidences]))
    
    all_algorithms <- c("COSG_All", "FindAllMarkers_All", "GenesorteR_All", "presto_All")
    all_confidences <- c("Conf_COSG_All", "Conf_FindAllMarkers_All", "Conf_GenesorteR_All", "Conf_presto_All")
    all_preds <- unlist(row[all_algorithms])
    all_confs <- as.numeric(unlist(row[all_confidences]))
    
    # 处理NA值
    semi_preds[is.na(semi_preds)] <- "Unknown"
    semi_confs[is.na(semi_confs)] <- 0
    all_preds[is.na(all_preds)] <- "Unknown"
    all_confs[is.na(all_confs)] <- 0
    
    # 获取半监督和无监督的最佳预测
    best_semi <- get_best_prediction(semi_preds, semi_confs)
    best_all <- get_best_prediction(all_preds, all_confs)
    
    # 如果半监督和无监督的最佳预测都存在且相同，则使用半监督规则
    if (!is.na(best_semi) && !is.na(best_all) && best_semi == best_all) {
      predicted_labels[i] <- process_semi_supervised(row)
      next
    }
    
    # 如果半监督和无监督的最佳预测不存在或不同，则按照以下规则
    # 规则3.1：所有置信度为0
    if (max(semi_confs) == 0 && max(all_confs) == 0) {
      predicted_labels[i] <- "Low quality"
      next
    }
    
    # 规则3.2：两种模式都有置信度大于0
    if (max(semi_confs) > 0 && max(all_confs) > 0) {
      # 如果半监督和无监督都有最佳预测，但不同
      if (!is.na(best_semi) && !is.na(best_all)) {
        # 特殊情况处理
        if (is_special_case(best_semi, best_all)) {
          # 获取最佳预测的置信度
          best_semi_conf <- get_best_confidence(semi_preds, semi_confs, best_semi)
          best_all_conf <- get_best_confidence(all_preds, all_confs, best_all)
          if (best_semi_conf > best_all_conf) {
            predicted_labels[i] <- best_semi
          } else if (best_all_conf > best_semi_conf) {
            predicted_labels[i] <- best_all
          } else {
            predicted_labels[i] <- "Undetermined"
          }
          next
        }
        
        # 非特殊情况，检查是否同一大类
        pan_semi <- get_pan_category(best_semi)
        pan_all <- get_pan_category(best_all)
        
        # 如果属于不同大类，标记为Low quality
        if (!is.na(pan_semi) && !is.na(pan_all) && pan_semi != pan_all) {
          predicted_labels[i] <- "Low quality"
          next
        }
        
        # 如果属于相同大类，按照半监督规则处理
        predicted_labels[i] <- process_semi_supervised(row)
        next
      } else {
        # 如果半监督或无监督中有一个没有最佳预测，则按照半监督规则处理（因为半监督优先级高）
        predicted_labels[i] <- process_semi_supervised(row)
        next
      }
    }
    
    # 规则3.3：只有半监督有置信度大于0
    if (max(semi_confs) > 0 && max(all_confs) == 0) {
      predicted_labels[i] <- process_semi_supervised(row)
      next
    }
    
    # 规则3.4：只有无监督有置信度大于0
    if (max(semi_confs) == 0 && max(all_confs) > 0) {
      predicted_labels[i] <- process_unsupervised(row)
      next
    }
    
    predicted_labels[i] <- "Undetermined" # 默认情况
  }
  
  # 关闭进度条
  close(pb)
  
  # 添加预测标签列到数据框
  df$predicted.label.finest <- predicted_labels
  df <- cbind(predicted.label.finest = df$predicted.label.finest, df)
  df$predicted.label.finest <- NULL
  
  # 显示完成信息和统计
  message("\nAnnotation completed successfully!")
  message("Summary of predicted labels:")
  label_table <- table(predicted_labels)
  for (label in names(label_table)) {
    message(paste0("  ", label, ": ", label_table[label], " cells"))
  }
  
  return(df)
}

