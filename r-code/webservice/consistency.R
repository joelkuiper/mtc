consistency <- function(params)  {
    d <- params$network
    d <- d$data
    d <- do.call(rbind, lapply(d, as.data.frame))

    # Remove 1-arm studies
    sel <- sapply(d$study, function(study) {
        sum(d$study == study) > 1
    })
    d <- d[sel, ]

    if(!is.null(params$network$treatments)) {
        treatments <- lapply(params$network$treatments, unlist)
        network <- mtc.network(data=d, treatments=treatments, description=unlist(params$network$description))
    } else {
        network <- mtc.network(data=d, description=unlist(params$network$description))
    }

    factor <- if(is.null(params$factor)) 2.5 else params$factor
    n.chain <- if(is.null(params$n_chain)) 4  else params$n_chain
    n.adapt <- if(is.null(params$n_adapt)) 5000  else params$n_adapt
    n.iter <- if(is.null(params$n_iter)) 20000  else params$n_iter
    thin <- if(is.null(params$thin)) 1  else params$thin

    model <- mtc.model(network, "consistency",  factor, n.chain)
    run <- mtc.run(model, n.adapt=n.adapt, n.iter=n.iter, thin=thin)

    rewrite.quantiles <- function(quantiles) {
      lapply(rownames(quantiles), function(param) {
          data <- as.list(quantiles[param, ])
          data$parameter <- param
          data
      })
    }

    quantiles <- rewrite.quantiles(summary(run)$quantiles)
    psrf <- gelman.diag(run)$psrf
    rank.prob <- rank.probability(run)

    # comprehensive list of relative effects
    calc.relative.effects <- function(result) {
        model <- result$model
        ts <- model$network$treatments$id
        n <- length(ts)
        idx <- rep(TRUE, times=n*n)
        f <- function(i) { (i - 1) * n + i }
        idx[f(1:n)] <- FALSE
        t1 <- rep(ts, each=length(ts))[idx]
        t2 <- rep(ts, times=length(ts))[idx]
        relative.effects <- summary(relative.effect(result, t1=t1, t2=t2, preserve.extra=F))$quantiles
        if (ll.call("scale.log", model)) {
            relative.effects <- exp(relative.effects)
        }
        relative.effects
    }
    relative.effects <- rewrite.quantiles(calc.relative.effects(run))

    save.plot(function() forest(run), "forest")
    save.plot(function() forest(run), "forest_vec", "svg")
    save.plot(function() plot(model), "model")
    save.plot(function() plot(model), "model_vec", "svg")
    ranks <- function() barplot(t(rank.prob),
                                col=rainbow(dim(rank.prob)[[1]]),
                                beside=T,
                                legend.text=paste("Rank", rep(1:dim(rank.prob)[[1]])))

    save.plot(ranks, "ranks")
    save.plot(ranks, "ranks_vec", "svg")

    results <- list("treatments" = network$treatments,
                    "description" = network$description,
                    "quantiles" = quantiles,
                    "psrf" = psrf,
                    "ranks" = rank.prob,
                    "relative_effects" = relative.effects)
    results
}
