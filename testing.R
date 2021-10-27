mammals <- read.table(
  "https://www.math.ntnu.no/~jarlet/statmod/mammals.dat",
  header=T)

mod0 <- lm(log(brain) ~ log(body), data = mammals)

is.human = ifelse(mammals$species == "Human", 1, 0)
mammals$is.human = as.factor(is.human)

mod1 <- lm(log(brain) ~ log(body) + is.human, data = mammals)

# Wald test linear:
C <- matrix(c(0, 1, 0), nrow = 1)
d <-  as.vector(3/4)
r <- 1
p <- 3
n <-  nrow(mammals)
beta.hat <- mod1$coefficients
s2 <- deviance(mod1)/(n-p)
X <- model.matrix(mod1)
XtX.inv <- solve(t(X) %*% X)

w <-  t((C %*% beta.hat - d)) %*% solve(s2*C %*% XtX.inv %*% t(C)) %*% (C %*% beta - d)
p.val <-pchisq(w, r, lower.tail = FALSE)
p.val.F <- pf(w, r, n - p, lower.tail = FALSE)
p.val.F
p.val

# LRT-test linear
mod1.offset <- lm(log(brain) ~ is.human, offset = 3/4*log(body), data = mammals)
A <- logLik(mod1.offset)
B <- logLik(mod1)

X.stat <- -2 * (as.numeric(A)-as.numeric(B))
p.val <- pchisq(X.stat, df = 1, lower.tail = FALSE)
p.val

# Wald test gamma 
mod.gamma <- glm(brain ~ log(body) + is.human, family = Gamma(link = "log"), data = mammals)

beta.hat <- as.vector(mod.gamma$coefficients)
w <- t(C %*% beta.hat - d) %*% solve(C %*% vcov(mod.gamma) %*% t(C)) %*% (C %*% beta.hat - d)

p.val <- pchisq(w, r, lower.tail = FALSE)
p.val

# LRT test gamma
mod.gamma.offset <- glm(brain ~ 1 + is.human, family = Gamma(link = "log"), offset = 3/4*log(body), data = mammals)
A <- logLik(mod.gamma.offset)
B <- logLik(mod.gamma)

X.stat <- -2 * (as.numeric(A)-as.numeric(B))
p.val <- pchisq(X.stat, df = 1, lower.tail = FALSE)
p.val

