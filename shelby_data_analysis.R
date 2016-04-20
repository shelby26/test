# look at sample zero-inflated lognormal distribution
install.packages("bbmle")
library(bbmle)

install.packages("fishMod")
library(fishMod)
mix
fm.dln <- deltaLN(ln.form=value~variable, binary.form=~variable, data=mix)
plot( fm.dln$fitted, fm.dln$residuals[,"quantile"], pch=20, main="Delta Log-Normal quantile residuals")
fm.dln$coef

exp(mean(log(c(1327.4, 5781.0, 2213.1, 3140.5, 526.0, 8629.8))))

#Random poisson distribution
data.pois = data.frame(Trt = c(rep("A", n), rep("B", n)), 
                       Response = c(rpois(n, mean.A), rpois(n, mean.B))
)
hist(data.pois)

library(MASS)
data.nb = data.frame(Trt = c(rep("A", n), rep("B", n)), 
                     Response=c(rnegbin(n, mean.A, 5), rnegbin(n, mean.B, 5))
)
hist(data.nb)







alone <- rbind(data.frame(variable="A", value=c(1545.3, 1303.2, 1479.1, 552.1, 2013.7, 24945.9, 8298.5, 2582.3, 5984.1, 25941.8, 2951.2, 18071.7)),
               data.frame(variable="B", value=c(704.7, 693.4, 11749.0, 5308.8, 1066.6, 4385.3, 4830.6)))
alone
hist(alone$value)

# Try regular Poisson distribution for mixture data
mod <- glm(round(value) ~ variable, family=poisson, data=alone)
summary(mod)
nd = data.frame(variable = c("A","B"))
nd
cbind(nd, 
      Mean = predict(mod, newdata=nd, type="response"), 
      SE = predict(mod, newdata=nd, type="response", se.fit=T)$se.fit
)

1 - pchisq(summary(mod)$deviance, 
           summary(mod)$df.residual
)   # Poisson does not fit the data


# Try negative binomial
model.nb = glm.nb(round(value) ~ variable, data = alone)
summary(model.nb)

1 - pchisq(summary(model.nb)$deviance,
           summary(model.nb)$df.residual
) # Negative binomial fits the data

cbind(nd, 
      Mean = predict(model.nb, newdata = nd, type = "response"), 
      SE = predict(model.nb, newdata = nd, type="response", se.fit = T)$se.fit
)




mix <- data.frame(A=c(1327.4, 5781.0, 0.0, 0.0, 2213.1, 0.0, 3140.5, 526.0, 8629.8, 0.0),
                  B=c(2301.4, 3013.0, 305.5, 273.5, 756.8, 3793.1, 7161.4, 1524.1, 3758.4, 1285.3))
mix <- melt(mix)
mix
hist(mix$value)

# Try regular Poisson distribution for mixture data
mod <- glm(round(value) ~ variable, family=poisson, data=mix)
summary(mod)
nd = data.frame(variable = c("A","B"))
nd
cbind(nd, 
      Mean = predict(mod, newdata=nd, type="response"), 
      SE = predict(mod, newdata=nd, type="response", se.fit=T)$se.fit
)

1 - pchisq(summary(mod)$deviance, 
           summary(mod)$df.residual
)   # Poisson does not fit the data


# Try ZERO-INFLATED POISSON DISTRIBUTION for mixture data
library(VGAM)

library(pscl)
model.zip = zeroinfl(round(value) ~ variable|variable, dist="poisson", data = mix)
?zeroinfl
summary(model.zip)

cbind(nd, 
      Count = predict(model.zip, newdata = nd, type = "count"),
      Zero = predict(model.zip, newdata = nd, type = "zero")
)

10^mean(log10(mix[which(mix$variable=="B"), "value"]))
mean(mix[which(mix$variable=="A"), "value"])
mean(c(1327.4, 5781.0, 2213.1, 3140.5, 526.0, 8629.8))

# Zero-inflated negative binomial
model.zip.3 = zeroinfl(round(value) ~ variable|variable, data = mix, dist = "negbin")
summary(model.zip.3)
# Theta parameter is not significant --> use zero-inflated poisson




