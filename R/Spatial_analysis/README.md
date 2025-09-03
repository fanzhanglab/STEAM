# STEAM Spatial Anchor Analysis - Complete Organization Guide

### Core Purpose
The method identifies and corrects misclassified cells by using high-confidence cells as "spatial anchors" to guide corrections of nearby uncertain cells, without requiring ground truth validation data.

### Scientific Logic

#### The Spatial Anchor Concept
- **Spatial Anchors**: High-confidence cells that serve as reference points for nearby corrections
- **Biological Basis**: Cells in close spatial proximity often belong to the same cell type due to tissue organization
- **Expression Validation**: Gene expression similarity provides additional evidence for corrections

#### Key Biological Principles
1. **Spatial Coherence**: Adjacent cells in tissue sections typically share cell type identity
2. **Expression Consistency**: Cells of the same type exhibit similar gene expression patterns  
3. **Confidence Weighting**: High-confidence predictions are more likely to be correct
4. **Consensus Voting**: Multiple spatial anchors reduce correction errors through voting

### Analysis Workflow

```
1. Compute Baseline Metrics
   - (Spatial coherence, expression coherence, model confidence)
   
2. Identify Uncertain Cells  
   - (Low confidence + spatial inconsistency + expression mismatch)
   
3. Find Spatial Anchors
   - (High-confidence cells in neighborhoods)
   
4. Generate Corrections via Consensus
   - (Anchor voting with spatial/expression validation)
   
5. Apply Conservative Corrections
   - (Only high-confidence improvements)
   
6. Iterate Until Convergence
   - (Repeat with updated predictions)
   
7. Track & Validate Improvements
```

### Mathematical Framework

#### Confidence Scoring
- **Spatial Confidence**: Agreement with k-nearest spatial neighbors
- **Expression Confidence**: Correlation with same-type cells  
- **Model Confidence**: Original classifier probability scores
- **Combined Score**: Weighted combination of all confidence measures

#### Correction Criteria
- **Consensus Threshold**: Minimum agreement among spatial anchors (typically >70%)
- **Improvement Threshold**: Minimum spatial/expression improvement required
- **Conservative Application**: Only corrections with high confidence are applied

### Key Advantages

1. **No Ground Truth Required**: Works without validation labels
2. **Conservative Approach**: Minimizes false corrections through multiple validation steps  
3. **Iterative Improvement**: Progressively refines predictions over multiple rounds
4. **Biologically Informed**: Leverages spatial tissue organization principles
5. **Comprehensive Validation**: Uses spatial, expression, and model confidence

---

## File Structure Overview

| File | Functions | Sections | Purpose |
|------|-----------|----------|---------|
| `iterativeHelpers.R` | 21 | 5 | Core helper functions for no-ground-truth iterative analysis |
| `STEAM_anchor.R` | 47 | - | Main STEAM anchor analysis implementation |
| `additional_helpers.R` | 21 | 4 | Advanced helper functions for enhanced iterative analysis |
| `plots.R` | 24 | - | Comprehensive visualization functions |
dashboard.R:             9 functions (Interactive dashboards)
foldSummary.R:           2 functions (Summary extraction)
accuracy_plots.R:        2 functions (Accuracy visualization)
----------------------------------------
TOTAL:                 126 functions
```

## 🔗 How the Components Work Together

### 🔄 Analysis Pipeline Flow

```
STEAM_anchor.R (Main Interface)
        ↓
    runIterative() 
        ↓
iterativeHelpers.R (Core Processing)
        ↓
1. baselineMetrics() → Establish validation baselines
2. identifyUncertain() → Find cells needing correction  
3. generateCorrections() → Create correction candidates
4. holdoutValidation() → Validate correction quality
5. applyCorrections() → Apply validated corrections
        ↓
additional_helpers.R (Advanced Processing)
        ↓
6. processEnhancedFoldClean() → Process individual folds
7. computeEnhancedConfidenceSimple() → Enhanced confidence
8. runIterationAnalysis() → Comprehensive iteration analysis
        ↓
Visualization & Summary
        ↓
9. plots.R → Spatial visualization & progress tracking
10. dashboard.R → Interactive exploration dashboards  
11. accuracy_plots.R → Performance visualization
12. foldSummary.R → Results summarization

```

### Component Interactions

#### **Core Analysis Loop**
1. **STEAM_anchor.R** provides the main `runIterative()` function
2. **iterativeHelpers.R** handles the core correction logic
3. **additional_helpers.R** processes folds and advanced analysis
4. Results flow back to update the STEAM object

#### **Data Flow**
```
STEAM Object → Spatial Coordinates + Expression Data
     ↓
Baseline Metrics (spatial coherence, expression coherence)
     ↓  
Uncertain Cell Identification (confidence + spatial + expression criteria)
     ↓
Spatial Anchor Analysis (high-confidence neighbors)
     ↓
Correction Generation (consensus voting among anchors)
     ↓
Validation & Application (holdout validation + conservative application)
     ↓
Updated STEAM Object → Visualization & Summary
```
