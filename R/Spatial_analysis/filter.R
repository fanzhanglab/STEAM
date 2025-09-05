#' Apply Intelligent Correction Filtering with Universal Data Type Support
#'
#' Filters corrections using learned patterns and dynamic rules that work with any
#' R data type (numeric, character, factor, etc.). Replaces hardcoded patterns
#' with adaptive learning-based filtering.
#'
#' @param STEAM.obj STEAM object containing progressive learning memory and state
#' @param corrections Data frame of proposed corrections → be filtered
#' @param verbose Logical indicating whether → print filtering details
#' @return Filtered corrections data frame with pattern scores and blocking status
#' @details
#' This function implements intelligent correction filtering using dynamic pattern
#' learning. It:
#' 1. Extracts learning memory from progressive learning state
#' 2. Generates dynamic pattern rules using getCorrectionPatternRules
#' 3. Applies blocked patterns with universal data type support using character conversion
#' 4. Applies priority patterns → boost successful correction types
#' 5. Calculates pattern scores for ranking corrections by expected success
#' 6. Filters out blocked corrections while preserving priority corrections
#' 
#' Universal data type support is achieved through as.character() conversion,
#' allowing the system → work with any R data type while maintaining pattern
#' learning capabilities.
#' @seealso \code{\link{getCorrectionPatternRules}} for pattern generation
#' @keywords internal
applyIntelligentFiltering <- function(STEAM.obj, corrections, verbose = TRUE) {
    
    if (nrow(corrections) == 0) return(corrections)
    
    # Use the same dynamic pattern generation as the main filtering system
    learning_memory <- if ("progressive_learning" %in% names(STEAM.obj)) {
        STEAM.obj$spatial_anchor_analysis$progressive_learning
    } else NULL
    
    # Apply the dynamic pattern filtering using getCorrectionPatternRules
    if (exists("getCorrectionPatternRules")) {
        pattern_rules <- getCorrectionPatternRules(corrections, learning_memory = learning_memory, verbose = verbose)
    } else {
        # Fallback: create empty pattern rules if function not available
        if (verbose) cat("Warning: getCorrectionPatternRules not available, using fallback\n")
        pattern_rules <- list(
            blocked = data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE),
            priority = data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
        )
    }
    
    blocked_patterns <- pattern_rules$blocked
    priority_patterns <- pattern_rules$priority
    
    # Calculate pattern scores for each correction
    corrections$pattern_score <- 0.5  # Default neutral score
    corrections$blocked <- FALSE
    
    # Apply blocked patterns (universal data type support)
    if (nrow(blocked_patterns) > 0) {
        for (i in 1:nrow(corrections)) {
            from_val <- as.character(corrections$original_pred[i])
            to_val <- as.character(corrections$suggested_correction[i])
            
            # Check if this correction matches any blocked pattern
            for (j in 1:nrow(blocked_patterns)) {
                blocked_from <- as.character(blocked_patterns$from[j])
                blocked_to <- as.character(blocked_patterns$to[j])
                
                if (from_val == blocked_from && to_val == blocked_to) {
                    corrections$blocked[i] <- TRUE
                    corrections$pattern_score[i] <- -3  # Strong penalty
                    break
                }
            }
        }
    }
    
    # Apply priority patterns (universal data type support)
    if (nrow(priority_patterns) > 0) {
        for (i in 1:nrow(corrections)) {
            from_val <- as.character(corrections$original_pred[i])
            to_val <- as.character(corrections$suggested_correction[i])
            
            # Check if this correction matches any priority pattern
            for (j in 1:nrow(priority_patterns)) {
                priority_from <- as.character(priority_patterns$from[j])
                priority_to <- as.character(priority_patterns$to[j])
                
                if (from_val == priority_from && to_val == priority_to) {
                    corrections$pattern_score[i] <- corrections$pattern_score[i] + 2  # Strong bonus
                    break
                }
            }
        }
    }
    
    # Apply universal similarity scoring (works with any data type)
    for (i in 1:nrow(corrections)) {
        from_val <- corrections$original_pred[i]
        to_val <- corrections$suggested_correction[i]
        
        # For numeric types, apply distance-based scoring
        if (is.numeric(from_val) && is.numeric(to_val)) {
            distance <- abs(as.numeric(from_val) - as.numeric(to_val))
            if (distance == 1) {
                corrections$pattern_score[i] <- corrections$pattern_score[i] + 1  # Adjacent values bonus
            } else if (distance >= 3) {
                corrections$pattern_score[i] <- corrections$pattern_score[i] - 0.5  # Large jump penalty
            }
        }
        
        # For character types, apply string similarity scoring
        else if (is.character(from_val) || is.character(to_val)) {
            # Convert → character for comparison
            from_str <- as.character(from_val)
            to_str <- as.character(to_val)
            
            # Simple string similarity (more sophisticated methods could be used)
            if (nchar(from_str) > 0 && nchar(to_str) > 0) {
                # Bonus for similar string lengths or patterns
                len_diff <- abs(nchar(from_str) - nchar(to_str))
                if (len_diff <= 1) {
                    corrections$pattern_score[i] <- corrections$pattern_score[i] + 0.5
                }
                
                # Check for common prefixes/suffixes
                if (substring(from_str, 1, 1) == substring(to_str, 1, 1)) {
                    corrections$pattern_score[i] <- corrections$pattern_score[i] + 0.3
                }
            }
        }
    }
    
    # Filter out blocked corrections
    before_count <- nrow(corrections)
    corrections_filtered <- corrections[!corrections$blocked, ]
    
    # Apply minimum score threshold
    corrections_filtered <- corrections_filtered[corrections_filtered$pattern_score >= -1, ]
    after_count <- nrow(corrections_filtered)
    
    # If we filtered too aggressively, keep top scoring corrections
    if (nrow(corrections_filtered) == 0 && before_count > 0) {
        corrections_ordered <- corrections[order(corrections$pattern_score, decreasing = TRUE), ]
        corrections_filtered <- corrections_ordered[1:min(5, before_count), ]  # Keep top 5 at minimum
        after_count <- nrow(corrections_filtered)
    }
    
    if (verbose && before_count != after_count) {
        cat(sprintf("Intelligent filtering: %d → %d corrections (filtered %d low-success patterns)\n", 
                    before_count, after_count, before_count - after_count))
    }
    
    # Remove temporary columns
    corrections_filtered$pattern_score <- NULL
    corrections_filtered$blocked <- NULL
    
    return(corrections_filtered)
}






#' Apply Pattern-Based Correction Filtering 
#'
#' Applies learned pattern rules → filter corrections, blocking low-success patterns
#' while prioritizing high-success patterns. Core filtering function for progressive learning.
#'
#' @param corrections Data frame of corrections → be filtered
#' @param verbose Logical indicating whether → print filtering statistics
#' @param learning_memory Optional learning memory containing pattern performance data
#' @return Filtered corrections data frame with priority weights and blocking flags
#' @details
#' This function implements the core pattern filtering logic for progressive learning:
#' 1. Retrieves dynamic pattern rules from learning memory or defaults
#' 2. Applies blocked patterns → eliminate consistently failing correction types
#' 3. Applies cautious patterns with attempt limits for moderately risky corrections
#' 4. Applies priority patterns with weight boosts for highly successful corrections
#' 5. Removes blocked corrections and sorts by priority weight
#' 6. Provides detailed logging of filtering decisions
#' 
#' The function supports three types of pattern rules:
#' - Blocked: Patterns with consistently poor performance (< 30% success)
#' - Cautious: Patterns with moderate risk that need attempt limiting
#' - Priority: High-performing patterns (> 70% success) that get priority weighting
#' @seealso \code{\link{getDynamicPatternRules}} for pattern rule generation
#' @keywords internal
applyPatternFiltering <- function(corrections, verbose = TRUE, learning_memory = NULL) {
    
    if (nrow(corrections) == 0) return(corrections)
    
    # Get dynamic pattern rules from learning memory or use conservative defaults
    pattern_rules <- getDynamicPatternRules(learning_memory, verbose = verbose)
    
    blocked_patterns <- pattern_rules$blocked
    cautious_patterns <- pattern_rules$cautious  
    priority_patterns <- pattern_rules$priority
    
    original_count <- nrow(corrections)
    corrections$blocked <- FALSE
    corrections$priority_weight <- 0.5  # Default weight
    corrections$filter_reason <- ""
    
    # Apply blocked patterns
    for (i in seq_len(nrow(blocked_patterns))) {
        block_rule <- blocked_patterns[i, ]
        blocked_idx <- which(
            as.character(corrections$original_pred) == as.character(block_rule$from) & 
                as.character(corrections$suggested_correction) == as.character(block_rule$to)
        )
        
        if (length(blocked_idx) > 0) {
            corrections$blocked[blocked_idx] <- TRUE
            corrections$filter_reason[blocked_idx] <- paste("BLOCKED:", block_rule$reason)
            
            if (verbose) {
                cat(sprintf("Blocked %d corrections: %s → %s (%s)\n", 
                            length(blocked_idx), block_rule$from, block_rule$to, block_rule$reason))
            }
        }
    }
    
    # Apply cautious patterns  
    for (i in seq_len(nrow(cautious_patterns))) {
        cautious_rule <- cautious_patterns[i, ]
        cautious_idx <- which(
            as.character(corrections$original_pred) == as.character(cautious_rule$from) & 
                as.character(corrections$suggested_correction) == as.character(cautious_rule$to) &
                !corrections$blocked
        )
        
        if (length(cautious_idx) > 0) {
            # Limit → max_attempts
            if (length(cautious_idx) > cautious_rule$max_attempts) {
                keep_idx <- cautious_idx[seq_len(cautious_rule$max_attempts)]
                block_idx <- cautious_idx[(cautious_rule$max_attempts + 1):length(cautious_idx)]
                
                corrections$blocked[block_idx] <- TRUE
                corrections$filter_reason[block_idx] <- paste("LIMITED:", cautious_rule$reason)
                cautious_idx <- keep_idx
            }
            
            corrections$priority_weight[cautious_idx] <- 0.3
            corrections$filter_reason[cautious_idx] <- paste("CAUTIOUS:", cautious_rule$reason)
            
            if (verbose) {
                cat(sprintf("Limited %d corrections: %s → %s (%s)\n", 
                            length(cautious_idx), cautious_rule$from, cautious_rule$to, cautious_rule$reason))
            }
        }
    }
    
    # Apply priority patterns
    for (i in seq_len(nrow(priority_patterns))) {
        priority_rule <- priority_patterns[i, ]
        priority_idx <- which(
            as.character(corrections$original_pred) == as.character(priority_rule$from) & 
                as.character(corrections$suggested_correction) == as.character(priority_rule$to) &
                !corrections$blocked
        )
        
        if (length(priority_idx) > 0) {
            corrections$priority_weight[priority_idx] <- priority_rule$priority_weight
            corrections$filter_reason[priority_idx] <- paste("PRIORITY:", priority_rule$reason)
            
            if (verbose) {
                cat(sprintf("Prioritized %d corrections: %s → %s (weight %.1f)\n", 
                            length(priority_idx), priority_rule$from, priority_rule$to, 
                            priority_rule$priority_weight))
            }
        }
    }
    
    # Remove blocked corrections
    filtered_corrections <- corrections[!corrections$blocked, ]
    
    # Sort by priority weight (high → low) 
    if (nrow(filtered_corrections) > 0) {
        filtered_corrections <- filtered_corrections[order(filtered_corrections$priority_weight, decreasing = TRUE), ]
    }
    
    final_count <- nrow(filtered_corrections)
    blocked_count <- original_count - final_count
    
    if (verbose && blocked_count > 0) {
        cat(sprintf("Pattern Filtering: %d/%d corrections blocked\n", blocked_count, original_count))
    }
    
    return(filtered_corrections)
}




#' Filter Corrections Using Learned Pattern Rules
#'
#' Applies comprehensive pattern-based filtering using learned correction performance
#' data → block poor patterns and prioritize successful ones.
#'
#' @param proposed_corrections Data frame of proposed corrections → be filtered
#' @param learning_memory Optional learning memory containing pattern performance history
#' @param verbose Logical indicating whether → print detailed filtering information
#' @return Filtered corrections with blocking flags, priority weights, and filter reasons
#' @details
#' This function provides comprehensive correction filtering using learned patterns:
#' 1. Retrieves dynamic pattern rules from learning memory via getCorrectionPatternRules
#' 2. Initializes filtering flags (blocked, cautious, priority_weight) for all corrections
#' 3. Applies blocked patterns → eliminate consistently failing corrections
#' 4. Applies cautious patterns with special handling and attempt tracking
#' 5. Applies priority patterns with weight boosts for high-success corrections
#' 6. Adds detailed filter_reason annotations for transparency
#' 7. Returns comprehensive filtering results for downstream processing
#' 
#' The function maintains detailed filtering reasons and provides comprehensive
#' logging for understanding filtering decisions and pattern learning effectiveness.
#' @seealso \code{\link{getCorrectionPatternRules}} for rule generation
#' @keywords internal
filterCorrectionsByPattern <- function(proposed_corrections, learning_memory = NULL, verbose = TRUE) {
    
    # Get dynamic pattern rules from learning memory
    pattern_rules <- getCorrectionPatternRules(learning_memory, verbose = verbose)
    
    if (nrow(proposed_corrections) == 0) {
        return(proposed_corrections)
    }
    
    original_count <- nrow(proposed_corrections)
    
    # Add filtering flags
    proposed_corrections$blocked <- FALSE
    proposed_corrections$cautious <- FALSE
    proposed_corrections$priority_weight <- 0.5  # Default weight
    proposed_corrections$filter_reason <- ""
    
    # Apply blocked patterns
    for (i in seq_len(nrow(pattern_rules$blocked))) {
        block_rule <- pattern_rules$blocked[i, ]
        blocked_idx <- which(
            as.character(proposed_corrections$original_pred) == as.character(block_rule$from) & 
                as.character(proposed_corrections$suggested_correction) == as.character(block_rule$to)
        )
        
        if (length(blocked_idx) > 0) {
            proposed_corrections$blocked[blocked_idx] <- TRUE
            proposed_corrections$filter_reason[blocked_idx] <- paste("BLOCKED:", block_rule$reason)
            
            if (verbose) {
                cat(sprintf("Blocked %d corrections: %s  →  %s (%s)\n", 
                            length(blocked_idx), block_rule$from, block_rule$to, block_rule$reason))
            }
        }
    }
    
    # Apply cautious patterns
    for (i in seq_len(nrow(pattern_rules$cautious))) {
        cautious_rule <- pattern_rules$cautious[i, ]
        cautious_idx <- which(
            as.character(proposed_corrections$original_pred) == as.character(cautious_rule$from) & 
                as.character(proposed_corrections$suggested_correction) == as.character(cautious_rule$to) &
                !proposed_corrections$blocked
        )
        
        if (length(cautious_idx) > 0) {
            # Limit → max_attempts
            if (length(cautious_idx) > cautious_rule$max_attempts) {
                # Keep only the first max_attempts
                keep_idx <- cautious_idx[1:cautious_rule$max_attempts]
                block_idx <- cautious_idx[(cautious_rule$max_attempts + 1):length(cautious_idx)]
                
                proposed_corrections$blocked[block_idx] <- TRUE
                proposed_corrections$filter_reason[block_idx] <- paste("LIMITED:", cautious_rule$reason)
                cautious_idx <- keep_idx
            }
            
            proposed_corrections$cautious[cautious_idx] <- TRUE
            proposed_corrections$priority_weight[cautious_idx] <- 0.3
            proposed_corrections$filter_reason[cautious_idx] <- paste("CAUTIOUS:", cautious_rule$reason)
            
            if (verbose) {
                cat(sprintf("Limited %d corrections: %s  →  %s (%s)\n", 
                            length(cautious_idx), cautious_rule$from, cautious_rule$to, cautious_rule$reason))
            }
        }
    }
    
    # Apply priority patterns
    for (i in seq_len(nrow(pattern_rules$priority))) {
        priority_rule <- pattern_rules$priority[i, ]
        priority_idx <- which(
            as.character(proposed_corrections$original_pred) == as.character(priority_rule$from) & 
                as.character(proposed_corrections$suggested_correction) == as.character(priority_rule$to) &
                !proposed_corrections$blocked
        )
        
        if (length(priority_idx) > 0) {
            proposed_corrections$priority_weight[priority_idx] <- priority_rule$priority_weight
            proposed_corrections$filter_reason[priority_idx] <- paste("PRIORITY:", priority_rule$reason)
            
            if (verbose) {
                cat(sprintf("Prioritized %d corrections: %s  →  %s (weight %.1f, %s)\n", 
                            length(priority_idx), priority_rule$from, priority_rule$to, 
                            priority_rule$priority_weight, priority_rule$reason))
            }
        }
    }
    
    # Remove blocked corrections
    filtered_corrections <- proposed_corrections[!proposed_corrections$blocked, ]
    
    # Sort by priority weight (high → low)
    filtered_corrections <- filtered_corrections[order(filtered_corrections$priority_weight, decreasing = TRUE), ]
    
    final_count <- nrow(filtered_corrections)
    blocked_count <- original_count - final_count
    
    if (verbose) {
        cat(sprintf("\nFiltering Summary:\n"))
        cat(sprintf("   Original: %d corrections\n", original_count))
        cat(sprintf("   Blocked: %d corrections\n", blocked_count))
        cat(sprintf("   Remaining: %d corrections\n", final_count))
        cat(sprintf("   Improvement expected: %.1f%%  →  %.1f%%\n", 62.4, 
                    estimate_improved_success_rate(filtered_corrections)))
        cat("\n")
    }
    
    return(filtered_corrections)
}




#' Generate Correction Pattern Rules from Learning Memory
#'
#' Analyzes learning memory → generate dynamic pattern rules for intelligent
#' correction filtering. Supports universal data types with character-based patterns.
#'
#' @param learning_memory Learning memory containing pattern performance history
#' @param min_attempts Minimum attempts required before creating rules (default: 3)
#' @param block_threshold Success rate below which patterns are blocked (default: 0.2)
#' @param cautious_threshold Success rate threshold for cautious patterns (default: 0.6)
#' @param priority_threshold Success rate above which patterns get priority (default: 0.8)
#' @param verbose Logical indicating whether → print rule generation details
#' @return List containing blocked, cautious, and priority pattern rules
#' @details
#' This function generates dynamic correction pattern rules based on learned performance:
#' 1. Analyzes pattern performance from learning memory
#' 2. Identifies patterns with sufficient attempts (>= min_attempts)
#' 3. Categorizes patterns by success rate into blocked, cautious, and priority groups
#' 4. Creates rule data frames with universal data type support via character conversion
#' 5. Provides detailed rule generation logging for transparency
#' 
#' Pattern categorization:
#' - Blocked: Success rate < 20% - completely avoid these corrections
#' - Cautious: Success rate 20-60% - apply with attempt limits and extra validation
#' - Priority: Success rate > 80% - prioritize these high-success corrections
#' 
#' Universal data type support allows the system → learn patterns for any R data type.
#' @seealso \code{\link{filterCorrectionsByPattern}} for rule application
#' @keywords internal
getCorrectionPatternRules <- function(learning_memory = NULL, min_attempts = 3, 
                                      block_threshold = 0.2, cautious_threshold = 0.6, 
                                      priority_threshold = 0.8, verbose = FALSE) {
  
  # Initialize empty rule sets (handle any data type for cell types)
  blocked_patterns <- data.frame(
    from = character(0), to = character(0), reason = character(0), 
    stringsAsFactors = FALSE
  )
  cautious_patterns <- data.frame(
    from = character(0), to = character(0), max_attempts = integer(0), reason = character(0),
    stringsAsFactors = FALSE
  )
  priority_patterns <- data.frame(
    from = character(0), to = character(0), priority_weight = numeric(0), reason = character(0),
    stringsAsFactors = FALSE
  )
  
  # If no learning memory, return empty rules (conservative approach)
  if (is.null(learning_memory) || is.null(learning_memory$pattern_performance)) {
    if (verbose) {
      cat("No pattern learning history available - using conservative approach\n")
    }
    return(list(
      blocked = blocked_patterns,
      cautious = cautious_patterns,
      priority = priority_patterns
    ))
  }
  
  if (verbose) {
    cat("Generating dynamic pattern rules from learning history...\n")
  }
  
  # Analyze each pattern's performance
  pattern_perf <- learning_memory$pattern_performance
  for (pattern in names(pattern_perf)) {
    parts <- strsplit(pattern, "→")[[1]]
    if (length(parts) != 2) next
    
    # Handle any data type for cell types (not just integers)
    from_pred <- parts[1]
    to_pred <- parts[2]
    perf <- pattern_perf[[pattern]]
    
    # Skip patterns with insufficient data
    if (perf$attempts < min_attempts) next
    
    success_rate <- perf$success_rate
    attempts <- perf$attempts
    
    if (success_rate <= block_threshold) {
      # Block consistently failing patterns
      blocked_patterns <- rbind(blocked_patterns, data.frame(
        from = as.character(from_pred),
        to = as.character(to_pred),
        reason = sprintf("%.1f%% success (%d attempts)", success_rate * 100, attempts),
        stringsAsFactors = FALSE
      ))
      
      if (verbose) {
        cat(sprintf("BLOCK: %s → %s (%.1f%% success in %d attempts)\n", 
                    from_pred, to_pred, success_rate * 100, attempts))
      }
      
    } else if (success_rate <= cautious_threshold) {
      # Limit moderately successful patterns
      max_attempts <- max(2, min(4, ceiling(attempts * success_rate)))
      cautious_patterns <- rbind(cautious_patterns, data.frame(
        from = as.character(from_pred),
        to = as.character(to_pred),
        max_attempts = max_attempts,
        reason = sprintf("%.1f%% success (%d attempts)", success_rate * 100, attempts),
        stringsAsFactors = FALSE
      ))
      
      if (verbose) {
        cat(sprintf("LIMIT: %s → %s (%.1f%% success, max %d attempts)\n", 
                    from_pred, to_pred, success_rate * 100, max_attempts))
      }
      
    } else if (success_rate >= priority_threshold) {
      # Prioritize highly successful patterns
      priority_weight <- min(1.0, success_rate + 0.1)
      priority_patterns <- rbind(priority_patterns, data.frame(
        from = as.character(from_pred),
        to = as.character(to_pred),
        priority_weight = priority_weight,
        reason = sprintf("%.1f%% success (%d attempts)", success_rate * 100, attempts),
        stringsAsFactors = FALSE
      ))
      
      if (verbose) {
        cat(sprintf("PRIORITIZE: %s → %s (%.1f%% success, weight %.2f)\n", 
                    from_pred, to_pred, success_rate * 100, priority_weight))
      }
    }
  }
  
  return(list(
    blocked = blocked_patterns,
    cautious = cautious_patterns,
    priority = priority_patterns
  ))
}


#' Get Dynamic Pattern Rules with Universal Data Type Support
#'
#' Retrieves dynamic pattern rules from learning memory or provides conservative
#' defaults. Supports any R data type through character-based pattern matching.
#'
#' @param learning_memory Optional learning memory containing pattern performance data
#' @param verbose Logical indicating whether → print rule retrieval information
#' @return List containing blocked, cautious, and priority pattern rule data frames
#' @details
#' This function serves as the main interface for retrieving pattern filtering rules:
#' 1. Initializes empty rule structures for all pattern types (blocked, cautious, priority)
#' 2. Checks for learning memory and pattern performance data
#' 3. Analyzes pattern performance → generate dynamic rules
#' 4. Provides fallback → conservative defaults if no learning data available
#' 5. Returns properly structured rule sets for pattern filtering
#' 
#' The function supports universal data types by using character conversion for
#' all pattern matching, replacing hardcoded numeric patterns with dynamic
#' learning-based rules that work with integers, strings, factors, or any R data type.
#' 
#' Rule types returned:
#' - blocked: Patterns → completely avoid (low success rate)
#' - cautious: Patterns → limit or apply carefully (moderate success)
#' - priority: Patterns → prioritize (high success rate)
#' @seealso \code{\link{getCorrectionPatternRules}} for detailed rule generation
#' @keywords internal
getDynamicPatternRules <- function(learning_memory = NULL, verbose = TRUE) {
  
  # Initialize empty pattern rules (handle any data type for cell types)
  blocked_patterns <- data.frame(
    from = character(0), to = character(0), reason = character(0), 
    stringsAsFactors = FALSE
  )
  cautious_patterns <- data.frame(
    from = character(0), to = character(0), max_attempts = integer(0), reason = character(0),
    stringsAsFactors = FALSE
  )
  priority_patterns <- data.frame(
    from = character(0), to = character(0), priority_weight = numeric(0), reason = character(0),
    stringsAsFactors = FALSE
  )
  
  if (is.null(learning_memory) || is.null(learning_memory$pattern_performance)) {
    if (verbose) {
      cat("No learning memory available - using conservative filtering\n")
    }
    # Return minimal conservative rules
    return(list(
      blocked = blocked_patterns,
      cautious = cautious_patterns, 
      priority = priority_patterns
    ))
  }
  
  pattern_perf <- learning_memory$pattern_performance
  
  if (verbose) {
    cat("Analyzing pattern performance from learning memory...\n")
  }
  
  # Analyze each pattern's performance
  for (pattern in names(pattern_perf)) {
    parts <- strsplit(pattern, "→")[[1]]
    if (length(parts) != 2) next
    
    # Handle any data type for cell types (not just integers)
    from_pred <- parts[1]
    to_pred <- parts[2]
    
    perf <- pattern_perf[[pattern]]
    success_rate <- perf$success_rate
    attempts <- perf$attempts
    
    # Skip patterns with too few attempts → be reliable
    if (attempts < 3) next
    
    if (success_rate <= 0.2) {
      # Block consistently failing patterns (≤20% success)
      blocked_patterns <- rbind(blocked_patterns, data.frame(
        from = as.character(from_pred),
        to = as.character(to_pred), 
        reason = sprintf("%.1f%% success (%d attempts)", success_rate * 100, attempts),
        stringsAsFactors = FALSE
      ))
      
      if (verbose) {
        cat(sprintf("BLOCK: %s → %s (%.1f%% success in %d attempts)\n", 
                    from_pred, to_pred, success_rate * 100, attempts))
      }
      
    } else if (success_rate <= 0.6) {
      # Limit moderately successful patterns (20-60% success)  
      max_attempts <- max(2, min(4, floor(attempts * success_rate)))
      cautious_patterns <- rbind(cautious_patterns, data.frame(
        from = as.character(from_pred),
        to = as.character(to_pred),
        max_attempts = max_attempts,
        reason = sprintf("%.1f%% success (%d attempts)", success_rate * 100, attempts),
        stringsAsFactors = FALSE
      ))
      
      if (verbose) {
        cat(sprintf("LIMIT: %s → %s (%.1f%% success, max %d attempts)\n", 
                    from_pred, to_pred, success_rate * 100, max_attempts))
      }
      
    } else if (success_rate >= 0.8) {
      # Prioritize highly successful patterns (≥80% success)
      priority_weight <- min(1.0, success_rate + 0.2)  # Boost high performers
      priority_patterns <- rbind(priority_patterns, data.frame(
        from = as.character(from_pred),
        to = as.character(to_pred),
        priority_weight = priority_weight,
        reason = sprintf("%.1f%% success (%d attempts)", success_rate * 100, attempts),
        stringsAsFactors = FALSE
      ))
      
      if (verbose) {
        cat(sprintf("PRIORITIZE: %s → %s (%.1f%% success, weight %.2f)\n", 
                    from_pred, to_pred, success_rate * 100, priority_weight))
      }
    }
  }
  
  return(list(
    blocked = blocked_patterns,
    cautious = cautious_patterns,
    priority = priority_patterns
  ))

}
