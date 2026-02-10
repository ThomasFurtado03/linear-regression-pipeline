####Assignment 1
####Thomas Furtado

#Packages----
library(ggplot2)
library(lmtest)
library(GGally)
library(car)
library(summarytools)


#Problem 1----
#i)----
X <- as.matrix(read.csv("model.matrix.csv"))

dim(X)
colnames(X)

#Build xtx
xtx <- t(X) %*% X

#a)
eigen_results <- eigen(xtx)

eigenvalues <- eigen_results$values
eigenvectors <- eigen_results$vectors

print(eigenvalues)
print(eigenvectors)

#b)

#The eigen values of xtx include one value which is very close to 0, -1.345885e-15. This means that the design matrix, X, is not full rank. Meaning that there is 1 column (out of 5) which is redundant.

colnames(X)
#The column "statetreated" seems to be linearly dependant on the x-intercept and on "stateuntreated". Therefore one predictor is redundant, resulting in multicollinearity.


#c)

test_eigenvector <- eigenvectors[, 1]
test_eigenvalue <- eigenvalues[1]

#normalizing eigenvector
norm_eigenvector <- test_eigenvector / sqrt(sum(test_eigenvector^2))

#Using Av = λv
#Here A = xtx

check_result <- xtx %*% norm_eigenvector
expected_result <- test_eigenvalue * norm_eigenvector

is_eigenvector <- all.equal(as.numeric(check_result), as.numeric(expected_result))
print(is_eigenvector)

#To verify an eigenvector, I checked the equation Av = λv, where A = xtx. Since Av was equal to λv, this confirms that this is an eigenvector of the matrix xtx.


#d)

df_X <- as.data.frame(X)
head(df_X)

#Lets add the response variable "rate"
data(Puromycin)
df_X$rate <- Puromycin$rate

#Lets fit the model, then evaluate it
fit <- lm(rate ~ ., data = df_X)
summary(fit)
coef(fit)
#Upon analysis, 2 coefficients could not be defined because of singularities. Also, X.Intercept AND statetreated have coefficients of NA. The model fit was attempted, but the coefficients could not be estimated.

# X.Intercept is shown as NA because there is redundancy in intercepts, as lm already incorperates an intercept, so there is a duplication, resulting in NA here.

#The column "statetreated" is, as recognized before, redundant and depends on other columns. This explains the NA observed here as well.

#The model was fitted, but NA values were returned for 2 coefficients due to singularities. This happened because the design matrix has redundant predictors. X.intercept duplicates interepts already introduced by lm(). And statetreated is linearly dependant on statuntreated. This means that the disng matrix is not full rank, so therfore all the coefficients cannot be estimated.

#Lets modify X (X_star) by eliminating redundant columns.
colnames(df_X)
#1 and 5 are the redundant ones

X_star <- df_X[, -c(1,5)]
colnames(X_star)

fit_star <- lm(rate ~ ., data = X_star)
summary(fit_star)
coef(fit_star)


#By modifying the design matrix, X, I removed the redundant columns and attempted to re-fit the model. The resulting model (fit_star) produced coefficients with no missing values, meaning it is much more fit, and a plausible model. The edited design matrix, X_star, is now full rank (no more collinearity.


#ii)----

#Lets examine the relationship.
data(Puromycin)

ggplot(Puromycin, aes(x = conc, y = rate)) +
  geom_point()

ggplot(Puromycin, aes(x = conc, y = rate, shape = state, color = state)) +
  geom_point(cex = 2)

#There is a positive relationship observed, as conc increases, so does rate.
#However, the relationship is not linear. So a straight-line model likely wont be appropriate.


#Lets anyways attempt to fit a basic model.
f0 <- lm(rate ~ conc + state, data = Puromycin)
summary(f0)
round(coef(f0), 3)

#Given the small p-value, concentration is a positive predictor or rate (p < 0.001). The residual standard error is large (about 26.7). The adjusted R^2 is about 0.69, therefore, the relationship may not be fully captured using the raw concentration scale.


#Lets try a log transformation
f1 <- lm(rate ~log(conc) + state, data = Puromycin)
summary(f1)
round(coef(f1), 3)

#The model fit improves quite a bit. Both log(conc) and state are significant predictors( p << 0.001) And the residual standard error drops down to about 11.4 (from 26.7). The adjusted R^2 also increases (0.95), informing me that log(conc) provides a much stronger linear relationship to rate, than conc on its original scale.


#visualize log
ggplot(Puromycin, aes(x = log(conc), y = rate, shape = state, color = state)) +
  geom_point(cex = 2)
#This now looks like a linear relationship!

#f1 assumes the same slope for both treatment groups, however the effect of conc may differ by group. Let's add an interaction term to investigate this.

#Lets try using an interaction term in the model, to see if the slope differs between the 2 state levels.
f2 <- lm(rate ~ log(conc) + state + log(conc):state, data = Puromycin)
summary(f2)
round(coef(f2), 3)

#Including the interaction term actually further imporves the model! The interaction coeffecient (log(conc):state) is statistically significant (p= 0.0027). This suggests that the relationship between conc and rate differs between treated and untreated groups!
#The residual standard error decreased once more (to 9.15), AND the adjusted R^2 increases to 0.963.

#Therefore, this interaction model seems to provide the best fit among the 3 tested models.



###########Equation for the f2 model
round(coef(f2), 3)

#rate = 209.194 + 37.110log(conc) - 44.606I -10.128[log(conc)I]

#Where, if state = treated, I = 0. And if state = untreated, I = 1



#This can also be written as 2 separate lines
#1) State = treated, I = 0 ------> rate = 209.194 + 37.110log(conc)

#2) State = untreated, I = 1 ----> rate = 209.194 - 44.606 + 37.110log(conc) - 10.128log(conc)
#                                  rate = 164.588 + 26.982log(conc)


#iii)----

#Lets examine the quality of fit using visualizations.

par(mar = c(5, 4, 4, 2) + 0.1)

#LINEARITY & HOMOSCEDASTICITY


#Residuals vs Fitted (linearity and variance)
plot(f2, which = 1,
     main = "Residuals vs Fitted\nLinearity & Homoscedasticity")
abline(h = 0, col = "red", lty = 2)

res_vs_fitted <- resid(f2)
fitted_values <- fitted(f2)
lines(lowess(fitted_values, res_vs_fitted), col = "blue", lwd = 2)

#Scale-Location plot (Homoscedasticity)
plot(f2, which = 3,
     main = "Scale-Location Plot\nConstant Variance Check")

#Residuals vs Predictor
residuals <- resid(f2)

plot(log(Puromycin$conc), residuals,
     main = "Residuals vs log(conc)\nHomoscedasticity Check",
     xlab = "log(conc)", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)
lines(lowess(log(Puromycin$conc), residuals), col = "blue", lwd = 2)


#NORMALITY

#Q-Q plot
plot(f2, which = 2,
     main = "Normal Q-Q Plot\nNormality Check")

#Histogram of residuals + normal curve
hist(residuals,
     main = "Histogram of Residuals",
     xlab = "Residuals",
     col = "lightblue",
     freq = FALSE)

x <- seq(min(residuals), max(residuals), length = 100)
y <- dnorm(x, mean = mean(residuals), sd = sd(residuals))
lines(x, y, col = "red", lwd = 2)

#Shapiro-Wilk test
shapiro_test <- shapiro.test(residuals)
cat("Shapiro-Wilk normality test:\n")
cat("W =", round(shapiro_test$statistic, 4),
    ", p-value =", format.pval(shapiro_test$p.value, digits = 4), "\n")

#INDEPENDANCE

#Residuals vs observation order
plot(residuals, type = "b",
     main = "Residuals vs Observation Order\nIndependence Check",
     xlab = "Observation Order", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)


#ACF plot
acf(residuals,
    main = "ACF of Residuals\nAutocorrelation Check")

#Durbin_Watson test
dw_test <- dwtest(f2)
cat("\nDurbin-Watson test:\n")
cat("DW =", round(dw_test$statistic, 4),
    ", p-value =", format.pval(dw_test$p.value, digits = 4), "\n")

#CONSTANT VARIANCE & INFLUENCE

#Breusch-Pagan test
bp_test <- bptest(f2)
cat("\nBreusch-Pagan test:\n")
cat("BP =", round(bp_test$statistic, 4),
    ", p-value =", format.pval(bp_test$p.value, digits = 4), "\n")

#Cook's distance plot
cooksdistance <- cooks.distance(f2)
plot(cooksdistance, type = "h",
     main = "Cook's Distance\nInfluential Points",
     ylab = "Cook's Distance")
abline(h = 4/length(cooksdistance), col = "red", lty = 2)

#The Shapiro-Wilk test for normality was non-significant (p = 0.38), indicating that the reisduals are normally distributed. 
#The Durbin-Watson test showed no evidence of autocorrelation (p = 0.42)
#The Breusch-Pagan test was non significant (p = 0.29), suggesting there is no stron heteroscedasticity.
#Lastly, the Cook's distance did not indicate any highly influential observations.


#iv)----

#Lets use my best model, f2, to predict the rate at conc = 0.15

#To do this for both states, let me create a small data frame
newdata <- data.frame(
  conc = 0.15,
  state = factor(c("treated", "untreated"), levels = levels(Puromycin$state))
)

predict(f2, newdata = newdata)


#At conc = 0.15 ppm, the model predicts a rate of about 138.8 for treated and 113.4 for untreated samples.
#Since the regression includes log(conc), treatment effect, and an interaction term, this means that both the intercept and slope depend on the state. The treatment modifies the relationship between conc and reaction rate.



#Problem 3----

help(swiss)
data(swiss)
str(swiss)
head(swiss)


#i)----

ggpairs(swiss)
#The pairwise matrix indicated strongly negative associations between Fertility and Education. And it shows positive associations with Catholic and Infant mortality. This info can be important to consider when devising the regression models.

#Lets start with full multiple regression
fit_full <- lm(Fertility ~ ., data = swiss)
summary(fit_full)

extractAIC(fit_full)[2]

#Lets now reduce using poor (very high) p-values
pvals <- summary(fit_full)$coeff[,4]
pvals

vnames <- c("Fertility", rownames(summary(fit_full)$coeff)[pvals < .1])
vnames <- vnames[vnames != "(Intercept)"]
vnames

fit_p <- lm(Fertility ~ ., data = swiss[, vnames])
summary(fit_p)
extractAIC(fit_p)[2]

#Lets compare the AIC of full model, to the AIC of reduced model.
extractAIC(fit_full)[2]
extractAIC(fit_p)[2]
#According to AIC, the reduced model is slightly better.

#Lets try automatic model selection
select_back <- step(fit_full, direction = c("backward"))
summary(select_back)
extractAIC(select_back)[2]


#Comparing the AIC values, the manual and automatic selection were exactly equal. So lets proceed from here.

#The factors that significantly influence fertility are agriculture, education, catholic, and infant mortality. All but one (Examination) proved to be contributing factors.

#Lets look at a few other factor's relationship
ggplot(swiss, aes(x = Education, y = Fertility)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE)
#Clear negative relationship

ggplot(swiss, aes(x = Catholic, y = Fertility)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE)


#Education shows a strong negative association with fertility, consistent with being the most influential predictor in the final model!
#Catholic shows a weaker but positive association, although the scatter may indicate that other factors also contribute.


#Let's assess multicollinearity for the final reduced model
vif(select_back)

#All VIF vales were below 2.2, indicating there is minimal multicollinearity among the remaining predictors.

#ii)----

ggplot(swiss, aes(x = Catholic, y = Fertility)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE)

#The scatter plot shows an slight upward trend, suggesting a positive association between Catholics and Fertility rates. However the relationship may not be very strong, as there is great variation among the different catholic levels. It has some relationship but Fertility is likely more influenced by other factors.

#Lets take a look at these 2 clusters of catholic levels, there is clearly a cluster at "low" numbers (0 - 20) and another at high numbers (80+). With a mysterious cluster at the 50 region.

#Lets perform the split by median, to examine these identified clusters.
swiss$catholic_group <- as.factor(
  ifelse(swiss$Catholic < median(swiss$Catholic), "low", "high")
)
#Quickly look at the groups
freq(swiss$catholic_group, headings = FALSE)
tapply(swiss$Fertility, swiss$catholic_group, mean)
#High catholic group has a slightly higher observed mean.

#Let's look at which provinces fall into each group.
rownames(swiss)[swiss$catholic_group == "low"]
rownames(swiss)[swiss$catholic_group == "high"]

fit_group <- lm(Fertility ~ catholic_group, data = swiss)
summary(fit_group)
#The difference in average fertility between low/high catholic groups is not significant (p = 0.154)


ggplot(swiss, aes(x = catholic_group, y = Fertility)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, height = 0)
#The high-catholic group has a higher median and mean fertility than the low-catholic group.
#The boxplots overlap, suggesting no dramatic differences.
#My p-value revealed that the average fertility rates between the 2 groups do not differ significantly at the 5% confidence level.



#iii)----

#Lets fit a regression model, and perform diagnostics!

m1 <- lm(Fertility ~ Examination, data = swiss)
summary(m1)
#(p = 9.45e-07 ). Examination is strongly significant.
#Interestingly, Examination is significant in a simple regression model, but is not kept in the multiple regression model, perhaps because it overlaps with other predictors.


#Lets now perform the diagnostics

#LINEARITY

# Residuals vs Fitted: check for curvature or patterns
plot(m1, which = 1)
abline(h = 0, col = "red", lty = 2)
lines(lowess(fitted(m1), resid(m1)), col = "blue", lwd = 2)


#Scatter looks even around 0, no strong curved pattern, the trend line is quite flat. Linearity seems satisfied.

#NORMALITY

#Normal Q–Q plot
plot(m1, which = 2)
#Most points follow the diagonal line, though there are deviations at the extremes (Glane, Rive Droite)

#Shapiro test for normality
shapiro.test(resid(m1))
#p-value = 0.6071, > 0.05. There is no evidence of non-normality.


#HOMOSCEDASTICITY

#Scale–Location plot
plot(m1, which = 3)
#The line is quite horizontal, the spread does not seems to increase substantially.

#Breusch–Pagan test
bptest(m1)
#p-value = 0.6192, > 0.05, no evidence of violation, the constant variance assumption is satisfied.

#INDEPENDENCE

#Durbin–Watson test for autocorrelation
dwtest(m1)
#p- value = 1.759e-08, << 0.05, The independence assumption is violated! The residuals show strong positive correlation. This may suggest that there is a regional factor at play regarding provinces, rather than independence among provinces!

#Autocorrelation function of residuals
acf(resid(m1))
#Early lags exceed the confidence bounds early on, agrees with positive correlation in residuals.


#INFLUENTIAL POINTS

# Cook's Distance
cooksdist <- cooks.distance(m1)
plot(cooksdist, type = "h",
     main = "Cook's Distance",
     ylab = "Cook's Distance")

abline(h = 4/length(cooksdist), col = "red", lty = 2)
#There are a few provinces which exceed the threshold line, this indicates the potential for influential points.

#The model is likely stable, although a few provinces may have a disproportionate influence.


#DIAGNOSTIC SUMMARY

#The plots and tests indiacte that the regression model Fertility ~ Examination satisfies linearity, normality, and contance homoscedasticity. 

#However, the Durbin-Watson test strongly rejected independence (p = 1.759e-08). This suggested positive correlation in residuals. The Cook's distance identified a few provinces with potentially influential impact, but most provinces had an acceptable influence.


#iv)----
#Extra findings

#1)
#Pairwise relationships/Correlation structure

ggpairs(swiss)

#Notable findings:
#Fertility has a strong negative correlation with Education (-0.664)
#Fertility has a strong negative correlation with Examination (-0.646)
#Fertility is positively correlated with Catholic and Infant Mortality, 0.464 and 0.417 respectively.

#The matrix also shows strong correlations among predictors! Education and Examination are highly correlated (0.698)!
#This explains why Examination was dropped in the multivariable model! It was decided that it overlaps with Education!



#2)
#Compare adj.R^2 across models

#Full model
summary(fit_full)$adj.r.squared

#Reduced model
summary(select_back)$adj.r.squared

#Simple model, Examination only
summary(m1)$adj.r.squared

#The adjusted R^2 values confirm that the reduced model provides almost identical power as the full model (Adj R^2 = 0.670714, and 0.670971 respecively).
#This indicates that Examination can be removed without any meaningful loss of fit.
#The Fertility ~ Examination showed a much lower Adj. R^2 (0.4042126), showing that its important to include multiple predictors to achieve an effective model.

#3)
#Lastly, lets look for multicollinearity in the full model.

vif(fit_full)
vif(select_back)

#In the full model, all VIF values are below 4, indicating no severe multicolinearity. However Examination showed a notable overlap with the other predictors (3.68).
#In the reduced model, all VIF values dropped. This means that multicollinearity is minimal, and the remaining predictors will provide information that is mostly independent.