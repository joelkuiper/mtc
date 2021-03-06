standardize.treatments <- function(treatments) {
    treatments$id <- as.factor(treatments$id)
    treatments$description <- as.character(treatments$description)
    rownames(treatments) <- as.character(treatments$id)
    treatments[order(treatments$id), ]
}

standardize.data <- function(data, treatment.levels) {
    data$study <- as.factor(data$study)
    data$treatment <- factor(as.character(data$treatment), levels=treatment.levels)
    data <- data[order(data$study, data$treatment), ]
    rownames(data) <- seq(1:nrow(data))
    data
}

mtc.network <- function(data, treatments=NULL, description="Network") {
    # standardize the data
    if (!is.data.frame(data)) { 
        data <- do.call(rbind, lapply(data, as.data.frame))
    }

    # standardize the treatments
    if (is.null(treatments)) {
        treatments <- unique(data$treatment)
    }
    if (is.list(treatments) && !is.data.frame(treatments)) { 
        treatments <- as.data.frame(do.call(rbind, treatments))
    }
    if (is.character(treatments) || is.factor(treatments)) {
        treatments <- data.frame(id=treatments, description=treatments)
    }
    treatments <- standardize.treatments(treatments)

    network <- list(
        description=description,
        treatments=treatments,
        data=standardize.data(data, levels(treatments$id))
    )

    mtc.network.validate(network)

    class(network) <- "mtc.network"
    network
}

read.mtc.network <- function(file) {
    doc <- XML::xmlInternalTreeParse(file)
    description <- unlist(XML::xpathApply(doc, "/network", XML::xmlGetAttr, "description"))
    type <- unlist(XML::xpathApply(doc, "/network", XML::xmlGetAttr, "type", "rate"))
    treatments <- XML::xpathApply(doc, "/network/treatments/treatment",
        function(node) {
            c(
                id = XML::xmlGetAttr(node, "id"),
                description = XML::xmlValue(node)
            )
        }
    )
    if (identical(type, "rate")) {
        data <- XML::xpathApply(doc, "/network/studies/study/measurement",
            function(node) {
                list(
                    study = XML::xmlGetAttr(XML::xmlParent(node), "id"),
                    treatment = XML::xmlGetAttr(node, "treatment"),
                    responders = as.numeric(XML::xmlGetAttr(node, "responders")),
                    sampleSize = as.numeric(XML::xmlGetAttr(node, "sample"))
                )
            }
        )
    } else if (identical(type, "continuous")) {
        data <- XML::xpathApply(doc, "/network/studies/study/measurement",
            function(node) {
                list(
                    study = XML::xmlGetAttr(XML::xmlParent(node), "id"),
                    treatment = XML::xmlGetAttr(node, "treatment"),
                    mean = as.numeric(XML::xmlGetAttr(node, "mean")),
                    std.dev = as.numeric(XML::xmlGetAttr(node, "standardDeviation")),
                    sampleSize = as.numeric(XML::xmlGetAttr(node, "sample"))
                )
            }
        )
    } else if (identical(type, "none")) {
        data <- XML::xpathApply(doc, "/network/studies/study/measurement",
            function(node) {
                list(
                    study = XML::xmlGetAttr(XML::xmlParent(node), "id"),
                    treatment = XML::xmlGetAttr(node, "treatment")
                )
            }
        )
    }
    mtc.network(data, treatments=treatments, description=description)
}

write.mtc.network <- function(network, file) {
    root <- XML::newXMLNode("network")
    XML::xmlAttrs(root)["description"] <- network$description
    type <- if ('responders' %in% colnames(network$data)) {
        'rate'
    } else if ('mean' %in% colnames(network$data)) {
        'continuous'
    } else {
        'none'
    }
    XML::xmlAttrs(root)["type"] <- type

    treatments <- XML::newXMLNode("treatments", parent = root)
    apply(network$treatments, 1, function(row) {
        node <- XML::newXMLNode("treatment", parent = treatments)
        XML::xmlAttrs(node)["id"] <- row['id']
        XML::xmlValue(node) <- row['description']
    })

    studies <- XML::newXMLNode("studies", parent = root)
    study <- sapply(levels(network$data$study), function(sid) {
        node <- XML::newXMLNode("study", parent = studies)
        XML::xmlAttrs(node)["id"] <- sid
        node
    })

    apply(network$data, 1, function(row) {
        node <- XML::newXMLNode('measurement', parent=study[[row['study']]])
        XML::xmlAttrs(node)['treatment'] <- row['treatment']
        if (identical(type, 'rate')) {
            XML::xmlAttrs(node)['responders'] <- row['responders']
            XML::xmlAttrs(node)['sample'] <- row['sampleSize']
        } else if (identical(type, 'continuous')) {
            XML::xmlAttrs(node)['mean'] <- row['mean']
            XML::xmlAttrs(node)['standardDeviation'] <- row['std.dev']
            XML::xmlAttrs(node)['sample'] <- row['sampleSize']
        }
        node
    })

    cat(XML::saveXML(root), file=file)
}

mtc.network.validate <- function(network) { 
    # Check that there is some data
    stopifnot(nrow(network$treatments) > 0)  
    stopifnot(nrow(network$data) > 0)

    # Check that the treatments are correctly cross-referenced and have valid names
    stopifnot(all(network$data$treatment %in% network$treatments$id))
    stopifnot(all(network$treatments$id %in% network$data$treatment))
    idok <- regexpr("^[A-Za-z0-9_]+$", network$treatment$id) != -1
    if(!all(idok)) {
        stop(paste('Treatment name "',
            network$treatment$id[which(!idok)], '" invalid.\n',
            ' Treatment names may only contain letters, digits, and underscore (_).'), sep='')
    }

    # Check that the data frame has a sensible combination of columns
    columns <- colnames(network$data)
    contColumns <- c('mean', 'std.dev', 'sampleSize')
    dichColumns <- c('responders', 'sampleSize')

    if (contColumns[1] %in% columns && dichColumns[1] %in% columns) {
        stop('Ambiguous whether data is continuous or dichotomous: both "mean" and "responders" present.')
    }

    if (contColumns[1] %in% columns && !all(contColumns %in% columns)) {
        stop(paste('Continuous data must contain columns:', paste(contColumns, collapse=', ')))
    }

    if (dichColumns[1] %in% columns && !all(dichColumns %in% columns)) {
        stop(paste('Dichotomous data must contain columns:', paste(dichColumns, collapse=', ')))
    }
}

as.treatment.factor <- function(x, network) {
    v <- network$treatments$id
    if (is.numeric(x)) {
        factor(x, levels=1:nlevels(v), labels=levels(v))
    } else if (is.factor(x) || is.character(x)) {
        x <- as.character(x)
        factor(x, levels=levels(v))
    }
}

mtc.study.design <- function(network, study) {
    data <- network$data
    sort(data$treatment[data$study == study])
}

coerce.factor <- function(x, prototype) {
    factor(x, levels=1:nlevels(prototype), labels=levels(prototype))
}

mtc.treatment.pairs <- function(treatments) {
    n <- length(treatments)
    t1 <- do.call(c, lapply(1:(n-1), function(i) { rep(treatments[i], n - i) }))
    t2 <- do.call(c, lapply(1:(n-1), function(i) { treatments[(i+1):n] }))
    data.frame(t1=coerce.factor(t1, treatments), t2=coerce.factor(t2, treatments))
}

# Get all comparisons with direct evidence from the data set.
# Returns a (sorted) data frame with two columns (t1 and t2).
mtc.comparisons <- function(network) {
    data <- network$data

    # Identify the unique "designs" (treatment combinations)
    design <- function(study) { mtc.study.design(network, study) }
    designs <- unique(lapply(levels(data$study), design))

    # Generate all pair-wise comparisons from each "design"
    comparisons <- do.call(rbind, lapply(designs, mtc.treatment.pairs))

    # Ensure the output comparisons are unique and always in the same order
    comparisons <- unique(comparisons)
    comparisons <- comparisons[order(comparisons$t1, comparisons$t2), ]
    row.names(comparisons) <- NULL
    comparisons$t1 <- as.treatment.factor(comparisons$t1, network)
    comparisons$t2 <- as.treatment.factor(comparisons$t2, network)
    comparisons
}

edges.create <- function(e, ...) {
    e <- t(matrix(c(e$t1, e$t2), ncol=2))
    edges(as.vector(e), ...)
}

graph.create <- function(v, e, ...) {
    g <- graph.empty()
    g <- g + vertex(levels(v))
    g <- g + edges.create(e, ...)
    g
}

mtc.network.graph <- function(network) {
    comparisons <- mtc.comparisons(network)
    treatments <- network$treatments$id
    graph.create(treatments, comparisons, arrow.mode=0)
}

## mtc.network class methods
print.mtc.network <- function(x, ...) {
    cat("MTC dataset: ", x$description, "\n", sep="")
    print(x$data)
}

summary.mtc.network <- function(object, ...) {
    studies <- levels(object$data[,1])
    m <- sapply(object$treatments[,1], function(treatment) {
        sapply(studies, function(study) { 
            any(object$data[,1] == study & object$data[,2] == treatment)
        })
    })
    colnames(m) <- object$treatments[,1]
    x <- as.factor(apply(m, 1, sum))
    levels(x) <- sapply(levels(x), function(y) { paste(y, "arm", sep="-") })
    list("Description"=paste("MTC dataset: ", object$description, sep=""),
         "Studies per treatment"=apply(m, 2, sum), 
         "Number of n-arm studies"=summary(x)) 
}

plot.mtc.network <- function(x, layout=igraph::layout.circle, ...) {
    igraph::plot.igraph(mtc.network.graph(x), layout=layout, ...)
}
