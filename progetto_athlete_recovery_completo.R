################################################################################
# PROGETTO R - Athlete Recovery
# Riassunto operativo dei 3 laboratori:
# 1) Regressione lineare, test, intervalli di confidenza, residui
# 2) Diagnostica: punti influenti, collinearità, normalità, Box-Cox
# 3) Selezione variabili, predizione, regressione logistica
################################################################################

##### 0. PACCHETTI #############################################################
# Installare solo se necessario:
# install.packages(c("car", "MASS", "leaps", "GGally", "corrplot", "caret", "pROC"))

library(car)
library(MASS)
library(leaps)
library(GGally)
library(corrplot)
library(caret)
library(pROC)

##### 1. IMPORTAZIONE E CONTROLLO DATI ########################################

# Mettere il file CSV nella working directory, oppure usare setwd("percorso/cartella")
# Controlla la cartella corrente:
getwd()

athlete_recovery <- read.csv("athlete_recovery_synthetic.csv")

# Visualizzazione iniziale
View(athlete_recovery)
head(athlete_recovery)
dim(athlete_recovery)
str(athlete_recovery)
summary(athlete_recovery)

# Controllo valori mancanti
sum(is.na(athlete_recovery))
colSums(is.na(athlete_recovery))

# Controllo nomi colonne
names(athlete_recovery)

# Conversione variabili categoriche
athlete_recovery$Athlete_ID    <- as.factor(athlete_recovery$Athlete_ID)
athlete_recovery$Gender        <- as.factor(athlete_recovery$Gender)
athlete_recovery$Sport_Type    <- as.factor(athlete_recovery$Sport_Type)
athlete_recovery$Training_Type <- as.factor(athlete_recovery$Training_Type)
athlete_recovery$Day_of_Week   <- as.factor(athlete_recovery$Day_of_Week)


# Conversione di Stress_Level in variabile numerica ordinale
# Low = 1, Medium = 2, High = 3
athlete_recovery$Stress_Level <- factor(
  athlete_recovery$Stress_Level,
  levels = c("Low", "Medium", "High"),
  ordered = TRUE
)

athlete_recovery$Stress_Level <- as.numeric(athlete_recovery$Stress_Level)
# Dataset pulito per i modelli: elimino righe con NA nelle variabili usate
athlete_clean <- na.omit(athlete_recovery)

################################################################################
##### 2. ANALISI ESPLORATIVA ###################################################
################################################################################

##### Variabili numeriche ###################################################

# Controlla struttura dataset
str(athlete_clean)

# Seleziona automaticamente SOLO le variabili numeriche
num_vars <- names(athlete_clean)[sapply(athlete_clean, is.numeric)]

num_vars

##### Scatterplot matrix ####################################################

pairs(athlete_clean[, num_vars], pch = 16)

# Scatterplot ridotto (più leggibile)
ggpairs(athlete_clean[, c("Recovery_Score",
                        "Training_Intensity",
                        "Training_Duration_Min",
                        "Sleep_Duration_Hours",
                        "Caffeine_Intake_mg",
                        "HRV_ms",
                        "Resting_Heart_Rate",
                        "Energy_Level")],
      pch = 16)

##### Correlation Matrix ####################################################

cor_matrix <- cor(athlete_clean[, num_vars],
                  use = "complete.obs")

round(cor_matrix, 2)

corrplot(cor_matrix,
         method = "color",
         type = "upper")

corrplot(cor_matrix,
         method = "number",
         type = "upper",
         tl.cex = 0.7)

##### Boxplot ###############################################################

boxplot(Recovery_Score ~ Gender,
        data = athlete_clean,
        main = "Recovery Score per genere",
        xlab = "Gender",
        ylab = "Recovery Score")

boxplot(Recovery_Score ~ Sport_Type,
        data = athlete_clean,
        main = "Recovery Score per sport",
        xlab = "Sport Type",
        ylab = "Recovery Score")

boxplot(Recovery_Score ~ Training_Type,
        data = athlete_clean,
        main = "Recovery Score per tipo di allenamento",
        xlab = "Training Type",
        ylab = "Recovery Score")

################################################################################
##### 3. REGRESSIONE LINEARE COMPLETA #########################################
################################################################################

# Modello completo: Recovery_Score come variabile risposta
mod_full <- lm(Recovery_Score ~ Age  + Training_Duration_Min + Training_Intensity +
                 Sleep_Duration_Hours + Caffeine_Intake_mg + Stress_Level +
                 Resting_Heart_Rate + HRV_ms + Mood_Score + Muscle_Soreness +
                 Energy_Level,
               data = athlete_clean)

summary(mod_full)

# Interpretazione:
# - Estimate: effetto stimato della variabile sul Recovery Score, a parità delle altre.
# - p-value: se < 0.05 la variabile è statisticamente significativa.
# - R-squared e Adjusted R-squared: quota di variabilità spiegata dal modello.
# - F-statistic: test globale del modello.

# Attributi utili del modello
mod_full$coefficients
mod_full$fitted.values
mod_full$residuals
vcov(mod_full)

################################################################################
##### 4. TEST F MANUALE E TEST T SUI COEFFICIENTI #############################
################################################################################


# Test F globale manuale: H0 = tutti i coefficienti, tranne intercetta, sono 0
SS_tot <- sum((athlete_clean$Recovery_Score - mean(athlete_clean$Recovery_Score))^2)
SS_res <- sum(mod_full$residuals^2)

p <- mod_full$rank
n <- nrow(model.frame(mod_full))

F_manual <- ((SS_tot - SS_res) / (p - 1)) / (SS_res / (n - p))
p_value_F <- 1 - pf(F_manual, p - 1, n - p)

F_manual
p_value_F

# Test t manuale per un coefficiente, esempio: HRV_ms
X <- model.matrix(mod_full)
sigma2 <- sum(mod_full$residuals^2) / (n - p)

coef_name <- "HRV_ms"
beta_hat <- coef(mod_full)[coef_name]
se_beta <- sqrt(sigma2 * diag(solve(t(X) %*% X))[coef_name])

T_manual <- abs(beta_hat / se_beta)
p_value_t <- 2 * (1 - pt(T_manual, n - p))

T_manual
p_value_t

################################################################################
##### 5. INTERVALLI DI CONFIDENZA #############################################
################################################################################

# Intervalli di confidenza automatici al 95%
confint(mod_full, level = 0.95)

# Intervallo manuale per HRV_ms
alpha <- 0.05
t_alpha2 <- qt(1 - alpha/2, n - p)
IC_HRV <- c(beta_hat - t_alpha2 * se_beta,
            beta_hat + t_alpha2 * se_beta)

################################################################################
##### 6. DIAGNOSTICA DEL MODELLO LINEARE ######################################
################################################################################

# Grafici diagnostici standard
par(mfrow = c(2, 2))
plot(mod_full)
par(mfrow = c(1, 1))

# Residui standardizzati
res_std <- mod_full$residuals / summary(mod_full)$sigma
plot(mod_full$fitted.values, res_std,
     xlab = "Fitted values", ylab = "Residui standardizzati",
     main = "Residui standardizzati vs fitted", pch = 16)
abline(h = c(-2, 2), lty = 2)

which(abs(res_std) > 2)

# Residui studentizzati
res_stud <- rstandard(mod_full)
plot(mod_full$fitted.values, res_stud,
     xlab = "Fitted values", ylab = "Residui studentizzati",
     main = "Residui studentizzati vs fitted", pch = 16)
abline(h = c(-2, 2), lty = 2)

which(abs(res_stud) > 2)

# QQ plot e Shapiro test: normalità dei residui
qqnorm(residuals(mod_full))
qqline(residuals(mod_full))
shapiro.test(sample(residuals(mod_full), min(5000, length(residuals(mod_full)))))
# Nota: con dataset grandi Shapiro è molto sensibile; guardare anche il QQ plot.

# Omoschedasticità: residui vs fitted
plot(mod_full$fitted.values, residuals(mod_full),
     xlab = "Fitted values", ylab = "Residui",
     main = "Residui vs Fitted", pch = 16)
abline(h = 0, lty = 2)

################################################################################
##### 7. PUNTI INFLUENTI: LEVERAGE E COOK #####################################
################################################################################

# Leverage points
lev <- hatvalues(mod_full)
threshold_lev <- 2 * 12 / n

plot(mod_full$fitted.values, lev,
     xlab = "Fitted values", ylab = "Leverage",
     main = "Plot dei leverage", pch = 16)
abline(h = threshold_lev, lty = 2)

id_lev <- which(lev > threshold_lev)
id_lev
length(id_lev)

points(mod_full$fitted.values[id_lev], lev[id_lev], pch = 16)

# Distanza di Cook
cook <- cooks.distance(mod_full)
threshold_cook <- 4 / (n - p)

plot(mod_full$fitted.values, cook,
     xlab = "Fitted values", ylab = "Cook distance",
     main = "Distanza di Cook", pch = 16)
abline(h = threshold_cook, lty = 2)

id_cook <- which(cook > threshold_cook)
id_cook
length(id_cook)

points(mod_full$fitted.values[id_cook], cook[id_cook], pch = 16)

# Influential plot complessivo
influencePlot(mod_full, main = "Influential Plot")

# Modello senza punti influenti secondo Cook, per confronto
athlete_no_cook <- athlete_clean[-id_cook, ]

mod_no_cook <- lm(formula(mod_full), data = athlete_no_cook)
summary(mod_no_cook)

# Confronto coefficienti tra modello completo e modello senza punti influenti
coef_comparison <- data.frame(
  Full_Model = coef(mod_full),
  No_Cook_Model = coef(mod_no_cook),
  Relative_Change = abs((coef(mod_no_cook) - coef(mod_full)) / coef(mod_full))
)
coef_comparison

################################################################################
##### 8. COLLINEARITÀ E VIF ###################################################
################################################################################
mod_full <- lm(Recovery_Score ~ Age  + Day + Training_Duration_Min + Training_Intensity +
                 Sleep_Duration_Hours + Caffeine_Intake_mg + Stress_Level +
                 Resting_Heart_Rate + HRV_ms + Mood_Score + Muscle_Soreness +
                 Energy_Level,
               data = athlete_clean)


summary(mod_full)
# VIF: valori sopra 5 indicano possibile collinearità; sopra 10 collinearità forte.
vif(mod_full)

# Se ci sono VIF alti, valutare di togliere variabili ridondanti.
# Esempio: se Day e Week risultano molto collineari, tenere solo una delle due.

mod_reduced_collinearity <- lm(Recovery_Score ~ Age + Training_Duration_Min + Training_Intensity +
                                 Sleep_Duration_Hours + Caffeine_Intake_mg + Stress_Level +
                                 Resting_Heart_Rate + HRV_ms + Mood_Score +
                                 Muscle_Soreness + Energy_Level,
                               data = athlete_clean)

summary(mod_reduced_collinearity)
vif(mod_reduced_collinearity)

################################################################################
##### 9. TRASFORMAZIONE BOX-COX ###############################################
################################################################################

# Box-Cox si può usare se la variabile risposta è positiva.
# Recovery_Score dovrebbe essere positiva, ma controlliamo.
min(athlete_clean$Recovery_Score)

# Se ci sono valori <= 0, aggiungere una costante positiva.
athlete_clean$Recovery_Positive <- athlete_clean$Recovery_Score
if(min(athlete_clean$Recovery_Positive) <= 0){
  athlete_clean$Recovery_Positive <- athlete_clean$Recovery_Positive - min(athlete_clean$Recovery_Positive) + 1
}

mod_boxcox_base <- lm(Recovery_Positive ~ Age + Day_of_Week + Training_Duration_Min +
                        Training_Intensity + Sleep_Duration_Hours + Caffeine_Intake_mg +
                        Stress_Level + Resting_Heart_Rate + HRV_ms + Mood_Score +
                        Muscle_Soreness + Energy_Level,
                      data = athlete_clean)

bc <- boxcox(mod_boxcox_base)
best_lambda <- bc$x[which.max(bc$y)]
best_lambda

# Se lambda è vicino a 1, la trasformazione non è utile.
# Se lambda è vicino a 0, provare log(y).

if(abs(best_lambda) < 0.15){
  mod_boxcox <- lm(log(Recovery_Positive) ~ Age + Day_of_Week + Training_Duration_Min +
                     Training_Intensity + Sleep_Duration_Hours + Caffeine_Intake_mg +
                     Stress_Level + Resting_Heart_Rate + HRV_ms + Mood_Score +
                     Muscle_Soreness + Energy_Level,
                   data = athlete_clean)
} else {
  athlete_clean$Recovery_BoxCox <- (athlete_clean$Recovery_Positive^best_lambda - 1) / best_lambda
  mod_boxcox <- lm(Recovery_BoxCox ~ Age + Day_of_Week + Training_Duration_Min +
                     Training_Intensity + Sleep_Duration_Hours + Caffeine_Intake_mg +
                     Stress_Level + Resting_Heart_Rate + HRV_ms + Mood_Score +
                     Muscle_Soreness + Energy_Level,
                   data = athlete_clean)
}

summary(mod_boxcox)
qqnorm(residuals(mod_boxcox))
qqline(residuals(mod_boxcox))

################################################################################
##### 10. SELEZIONE VARIABILI: AIC, BIC, R2 AGGIUSTATO #########################
################################################################################

# Modello di partenza per selezione: evito ID atleta e variabili troppo ridondanti.
mod_start <- lm(Recovery_Score ~ Age  + Training_Duration_Min + Training_Intensity +
                  Sleep_Duration_Hours + Caffeine_Intake_mg + Stress_Level +
                  Resting_Heart_Rate + HRV_ms + Mood_Score + Muscle_Soreness +
                  Energy_Level,
                data = athlete_clean)


# Backward selection con AIC
mod_aic <- step(mod_start, direction = "backward", trace = TRUE)
summary(mod_aic)
AIC(mod_start, mod_aic)
BIC(mod_start, mod_aic)

summary(mod_start)
BIC(mod_start)

# Backward selection con BIC
n_start <- nrow(model.frame(mod_start))
mod_bic <- step(mod_start, direction = "backward", k = log(n_start), trace = TRUE)
summary(mod_bic)
AIC(mod_bic)
BIC(mod_bic)

# Ricerca esaustiva con R2 aggiustato: solo variabili numeriche per semplicità
mod_numeric <- lm(Recovery_Score ~ Age + Training_Duration_Min + Training_Intensity +
                    Sleep_Duration_Hours + Caffeine_Intake_mg + Stress_Level +
                    Resting_Heart_Rate + HRV_ms + Mood_Score + Muscle_Soreness +
                    Energy_Level,
                  data = athlete_clean)

x <- model.matrix(mod_numeric)[, -1]
y <- athlete_clean$Recovery_Score

best_adjr2 <- leaps(x, y, method = "adjr2")
best_index <- which.max(best_adjr2$adjr2)
best_adjr2$which[best_index, ]
best_adjr2$adjr2[best_index]

################################################################################
##### 11. PREDIZIONE CON IL MODELLO LINEARE ###################################
################################################################################

# Divido in training e test set
set.seed(123)
train_index <- createDataPartition(athlete_clean$Recovery_Score, p = 0.8, list = FALSE)
train_data <- athlete_clean[train_index, ]
test_data  <- athlete_clean[-train_index, ]

mod_pred <- lm(formula(mod_aic), data = train_data)
summary(mod_pred)

# Predizioni puntuali sul test set
pred_test <- predict(mod_pred, newdata = test_data)

# Misure di errore
RMSE <- sqrt(mean((test_data$Recovery_Score - pred_test)^2))
MAE  <- mean(abs(test_data$Recovery_Score - pred_test))
R2_test <- cor(test_data$Recovery_Score, pred_test)^2

RMSE
MAE
R2_test

# Intervalli di confidenza per il valore medio E[Y|X]
pred_conf <- predict(mod_pred, newdata = test_data[1:10, ], interval = "confidence")
pred_conf

# Intervalli di predizione per nuove osservazioni singole
pred_interval <- predict(mod_pred, newdata = test_data[1:10, ], interval = "prediction")
pred_interval

# Grafico osservati vs predetti
plot(test_data$Recovery_Score, pred_test,
     xlab = "Recovery Score osservato", ylab = "Recovery Score predetto",
     main = "Osservati vs predetti", pch = 16)
abline(0, 1, lty = 2)


################################################################################
# 12. REGRESSIONE LOGISTICA - Creazione variabile binaria
################################################################################

athlete_clean$Recovery_High <- ifelse(
  athlete_clean$Recovery_Score >= 70,
  1,
  0
)

athlete_clean$Recovery_High <- as.factor(athlete_clean$Recovery_High)

table(athlete_clean$Recovery_High)
prop.table(table(athlete_clean$Recovery_High))

################################################################################
# 4. Regressione logistica completa
################################################################################

log_full <- glm(
  Recovery_High ~ Age +
    Training_Duration_Min +
    Training_Intensity +
    Sleep_Duration_Hours +
    Caffeine_Intake_mg +
    Stress_Level +
    Resting_Heart_Rate +
    HRV_ms +
    Mood_Score +
    Muscle_Soreness +
    Energy_Level,
  data = athlete_clean,
  family = binomial
)


summary(log_full)

################################################################################
# 5. Odds ratio e intervalli di confidenza
################################################################################

odds_ratio <- exp(coef(log_full))
odds_ratio

odds_ratio_ci <- exp(confint(log_full))
odds_ratio_ci

odds_table <- data.frame(
  Variable = names(odds_ratio),
  Odds_Ratio = as.numeric(odds_ratio),
  CI_lower = odds_ratio_ci[, 1],
  CI_upper = odds_ratio_ci[, 2]
)

odds_table

################################################################################
# 6. Grafico degli Odds Ratio
################################################################################

odds_table_plot <- odds_table %>%
  filter(Variable != "(Intercept)")

ggplot(
  odds_table_plot,
  aes(x = reorder(Variable, Odds_Ratio), y = Odds_Ratio)
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.2
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "Odds Ratio - Regressione logistica",
    x = "Variabile",
    y = "Odds Ratio"
  ) +
  theme_minimal()

################################################################################
# 7. Probabilità predette
################################################################################

athlete_clean$prob_recovery_high <- predict(
  log_full,
  type = "response"
)

summary(athlete_clean$prob_recovery_high)

################################################################################
# 8. Distribuzione delle probabilità predette
################################################################################

ggplot(
  athlete_clean,
  aes(x = prob_recovery_high, fill = Recovery_High)
) +
  geom_histogram(
    bins = 40,
    alpha = 0.6,
    position = "identity"
  ) +
  labs(
    title = "Distribuzione delle probabilità predette",
    x = "Probabilità stimata di recupero alto",
    y = "Frequenza",
    fill = "Recovery High"
  ) +
  theme_minimal()

################################################################################
# 9. Classificazione con soglia 0.5
################################################################################

athlete_clean$class_pred <- ifelse(
  athlete_clean$prob_recovery_high >= 0.5,
  1,
  0
)

athlete_clean$class_pred <- as.factor(athlete_clean$class_pred)

table(
  Predetto = athlete_clean$class_pred,
  Osservato = athlete_clean$Recovery_High
)

accuracy <- mean(
  athlete_clean$class_pred == athlete_clean$Recovery_High
)

accuracy



################################################################################
# 10. Curva ROC e AUC
################################################################################

roc_obj <- roc(
  athlete_clean$Recovery_High,
  athlete_clean$prob_recovery_high
)

plot(
  roc_obj,
  main = "Curva ROC - Regressione logistica",
  col = "blue",
  lwd = 2
)

auc_value <- auc(roc_obj)
auc_value
################################################################################
# VALIDAZIONE TRAIN/TEST - REGRESSIONE LOGISTICA
################################################################################

set.seed(123)

train_index_log <- createDataPartition(
  athlete_clean$Recovery_High,
  p = 0.8,
  list = FALSE
)

train_log <- athlete_clean[train_index_log, ]
test_log  <- athlete_clean[-train_index_log, ]

log_train <- glm(
  Recovery_High ~ Age +
    Training_Duration_Min +
    Training_Intensity +
    Sleep_Duration_Hours +
    Caffeine_Intake_mg +
    Stress_Level +
    Resting_Heart_Rate +
    HRV_ms +
    Mood_Score +
    Muscle_Soreness +
    Energy_Level,
  data = train_log,
  family = binomial
)

prob_test <- predict(
  log_train,
  newdata = test_log,
  type = "response"
)

class_test <- ifelse(prob_test >= 0.5, 1, 0)
class_test <- as.factor(class_test)

confusion_matrix_test <- table(
  Predetto = class_test,
  Osservato = test_log$Recovery_High
)

confusion_matrix_test

accuracy_test <- mean(class_test == test_log$Recovery_High)
accuracy_test

roc_test <- roc(
  test_log$Recovery_High,
  prob_test
)

plot(
  roc_test,
  main = "Curva ROC - Test set",
  col = "blue",
  lwd = 2
)

auc_test <- auc(roc_test)
auc_test

################################################################################
# 11. Grafico Recovery Score vs probabilità predetta
################################################################################

ggplot(
  athlete_clean,
  aes(x = Recovery_Score, y = prob_recovery_high)
) +
  geom_point(alpha = 0.4) +
  geom_vline(xintercept = 70, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(
    title = "Recovery Score e probabilità stimata di recupero alto",
    x = "Recovery Score",
    y = "Probabilità stimata"
  ) +
  theme_minimal()

##############################################################################

################################################################################
# FINE
################################################################################

