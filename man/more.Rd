% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/more.R
\name{more}
\alias{more}
\title{more: Multi-Omics Regulation}
\usage{
more(
  targetData,
  regulatoryData,
  associations = NULL,
  omicType = NULL,
  condition = NULL,
  clinic = NULL,
  clinicType = NULL,
  minVariation = 0,
  percNA = 0.2,
  scaleType = "auto",
  epsilon = 1e-05,
  interactions = TRUE,
  varSel = "EN",
  alfaEN = NULL,
  correlation = 0.7,
  groupingISGL = "MF_0.7",
  alfa = 0.05,
  vip = 0.8,
  method = "MLR",
  parallel = FALSE,
  seed = 123
)
}
\arguments{
\item{targetData}{Data frame containing data from a target omic with its features in rows and
samples (in different conditions) in columns. Row names must be the target omic features IDs; e.g. when Gene Expression is considered gene IDs.}

\item{regulatoryData}{List where each element corresponds to a different omic data type to be considered (miRNAs,
transcription factors, methylation, etc.). The names of the list will represent the omics, and each element in 
the list should be a data matrix with omic regulators in rows and samples in columns.}

\item{associations}{List where each element corresponds to a different omic data type (miRNAs,
transcription factors, methylation, etc.). The names of the list will represent the omics. Each element in 
the list should be a data frame with 2 columns (optionally 3), describing the potential interactions between target omic features
and regulators for that omic. First column must contain the target omic features (the ones of 
\code{\link{targetData}} object), second column must contain the regulators, and an optional third column can
be added to describe the type of interaction (e.g., for methylation, if a CpG site is located in
the promoter region of the gene, in the first exon, etc.). If the user lacks prior knowledge of the potential regulators, they can set the parameter to NULL. 
In this case, all regulators in \code{\link{regulatoryData}} will be treated as potential regulators for all target omic features. In this case, for computational efficiency, it is recommended to use PLS2 \code{\link{method}}.
Additionally, if the users have prior knowledge for certain omics and want to set other omics to NULL, they can do so.}

\item{omicType}{Vector which indicates the type of data of omics introduced in \code{\link{regulatoryData}}. The user should code as 0 numeric omics and as 1 categorical or binary omics. 
By default is set to NULL. In this case, the data type will be predicted automatically. However, the user must verify the prediction and manually input the vector if incorrect.}

\item{condition}{Data frame describing the condition or phenotype to which considered samples belong. Rows must be the samples (columns
in \code{\link{targetData}}) and columns must be the conditions or phenotypes to be included in the model (e.g., treatment, etc.).}

\item{clinic}{Data frame with all clinical variables to consider, with samples in rows and variables in columns.}

\item{clinicType}{Vector which indicates the type of data of variables introduced in \code{\link{clinic}}. The user should code as 0 numeric variables and as 1 categorical or binary variables. 
By default is set to NULL. In this case, the data type will be predicted automatically. However, the user must verify the prediction and manually input the vector if incorrect.}

\item{minVariation}{For numerical regulators, it specifies the minimum change required across conditions to retain the regulator in 
the regression models. In the case of binary regulators, if the proportion of the most common value is equal to or inferior this value, 
the regulator is considered to have low variation and will be excluded from the regression models. The user has the option to set a single 
value to apply the same filter to all omics, provide a vector of the same length as omics if they want to specify different levels for each omics, 
or use 'NA' when they want to apply a minimum variation filter but are uncertain about the threshold. By default, 0.}

\item{percNA}{Maximum percentage of missing values present in regulatoryData regulators and observations to use them to construct the models. Only used in PLS models.}

\item{scaleType}{Type of scaling to be applied. Four options:
\itemize{
\item auto : Applies the auto scaling. 
\item softBlock : Applies the pareto scaling. \deqn{\frac{X_k}{s_k \sqrt[4]{m_b}} }
\item hardBlock : Applies the block scaling. \deqn{ \frac{X_k}{s_k \sqrt{m_b}} }
\item none : It does not apply any type of scaling. Not recommended if the user does not apply their own scaling.
}
considering m_b the number of variables of the block. By default, auto.}

\item{epsilon}{Convergence threshold for coordinate descent algorithm in Elasticnet. Default value, 1e-5.}

\item{interactions}{If TRUE, the model includes interactions between regulators and condition variables. By default, TRUE.}

\item{varSel}{Type of variable selection method to apply, different options depending on MLR or PLS usage. Four options:
\itemize{
\item EN : Applies a Multiple Linear Regression (MLR) with ElasticNet regularization.
\item ISGL : Applies a Multiple Linear Regression (MLR) with Iterative Sparse Group Lasso (ISGL) regularization.
\item Jack : Applies Jack-Knife resampling technique for the calculation of the significance of the coefficients in Partial Least Squares (PLS) models.
\item Perm : Applies a resampling technique for the calculation of the significance of the coefficients in Partial Least Squares (PLS) models in which the response variable is permuted 100 times to obtain the distribution of the coefficients and compute then their associated p-value.
}
By default, EN.}

\item{alfaEN}{ElasticNet mixing parameter. There are three options:
\itemize{
\item NULL : The parameter is selected from a grid of values ranging from 0 to 1 with 0.1 increments. The chosen value optimizes the mean cross-validated error when optimizing the lambda values.
\item A number between 0 and 1 : ElasticNet is applied with this number being the combination between Ridge and Lasso penalization (elasticnet=0 is the ridge penalty, elasticnet=1 is the lasso penalty). 
\item A vector with the mixing parameters to try. The one that optimizes the mean cross-validated error when optimizing the lambda values will be used.
}
By default, NULL.}

\item{correlation}{Value to determine the presence of collinearity between two regulators when using the MLR \code{\link{method}}. By default, 0.7.}

\item{groupingISGL}{Type of approach to take for creating the groups for the Iterative Sparse Group Lasso regularization. There are two options:
\itemize{
\item MF_X: Takes the same approach as in the multicollinearity filter. Grouping the variables which correlate more than X. E.g: MF_0.7 groups variables in the same groups when they are correlated more than 0.7. By default, MF_0.7. 
\item PCA_X: Extracts as many PCA components to explain X percent of the variability of the data and assigns the variables to the component in which they had the highest loadings. 
}}

\item{alfa}{Significance level for variable selection in PLS1 and PLS2 \code{\link{method}}. By default, 0.05.}

\item{vip}{Value of VIP above which a variable can be considered significant in addition to the computed p-value in \code{\link{varSel}}. By default, 0.8.}

\item{method}{Model to be fitted. Three options:
\itemize{
\item MLR : Applies a Multiple Linear Regression (MLR).
\item PLS1 : Applies a Partial Least Squares (PLS) model, one for each of the features in the target omic of \code{\link{targetData}}.
\item PLS2 : Applies a PLS model to all features of the target omic simultaneously, only possible when \code{\link{associations}}= NULL.
}
By default, MLR.}

\item{parallel}{If FALSE, MORE will be run sequentially. If TRUE, MORE will be run using parallelization with as many cores as the available ones minus one so the system is not overcharged. If the user wants to specify how many cores they want to use, they can also provide the number of cores to use in this parameter. Parallelization is only implemented for MLR with EN variable selection and PLS methods.}

\item{seed}{Sets the seed to guaranty reproducibility of the results. By default, 123.}
}
\value{
List containing the following elements:
\itemize{
\item ResultsPerTargetF : List with as many elements as features of the target omic in \code{\link{targetData}}. For each feature, it includes information about the feature values, considered variables, estimated coefficients,
                   detailed information about all regulators, and regulators identified as relevant (in MLR scenario) or significant (in PLS scenarios).
\item GlobalSummary : List with information about the fitted models, including model metrics, information about regulators, features of the target omic without models, regulators, master regulators and hub target features.
\item Arguments : List containing all the arguments used to generate the models.                
}
}
\description{
\code{more} fits a MLR regression model (when the selected method is MLR) or a PLS model (when the selected method is PLS1 or PLS2) for all target features in the dataset to identify
the potential regulators that show a significant impact on the expression of the target omic features under specific experimental conditions.
}
\examples{

data(TestData)

#Omic type
omicType = c(1,0,0)
names(omicType) = names(TestData$data.omics)
SimMLR = more(targetData = TestData$GeneExpressionDE,
              associations = TestData$associations, 
              regulatoryData = TestData$data.omics,
              omicType = omicType,
              condition = TestData$edesign,
              scaleType = 'auto', varSel = 'EN', 
              epsilon = 0.00001, alfaEN = NULL,
              interactions = TRUE, minVariation = 0,  
              correlation = 0.7, method  ='MLR')

}
