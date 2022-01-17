library(GeneSummary)
library(tm)
library(org.Hs.eg.db)

## List of entrez ID
keggPathways <- org.Hs.egPATH2EG
mappedKeys <- mappedkeys(keggPathways)
keggList <- as.list(keggPathways[mappedKeys])
geneList <- keggList$`05160`
additionalRemove <- NA

## Using GeneSummary
tb <- loadGeneSummary()

## Already performed rda with high frequency words
allDocsPre <- VCorpus(VectorSource(tb$Gene_summary))
allDocsPre <- allDocsPre %>%
    tm_map(FUN=content_transformer(tolower)) %>% 
    tm_map(FUN=removeNumbers) %>%
    tm_map(removeWords, stopwords::stopwords("english", "stopwords-iso")) %>%
    tm_map(FUN=removePunctuation) %>%
    tm_map(FUN=stripWhitespace)
matAllPre <- as.matrix(TermDocumentMatrix(allDocsPre))
matAllPreSorted <- sort(rowSums(matAllPre), decreasing=TRUE)
allFreqGeneSummary <- data.frame(matAllPreSorted)
allFreqGeneSummary$word <- rownames(allFreqGeneSummary)
colnames(allFreqGeneSummary) <- c("freq","word")
# load("allFreqGeneSummary.rda") 

## Filter high frequency words
excludeFreq <- 5000
filterWords <- allFreqGeneSummary[allFreqGeneSummary$freq>excludeFreq,]$word
filterWords <- c(filterWords, "pmids", "geneid") # 'PMIDs' is excluded by default
fil <- tb %>% filter(Gene_ID %in% geneList)
   
## Make corpus for queried genes
docs <- VCorpus(VectorSource(fil$Gene_summary))
docs <- docs %>%
    tm_map(FUN=content_transformer(tolower)) %>% 
    tm_map(FUN=removeNumbers) %>%
    tm_map(removeWords, stopwords::stopwords("english", "stopwords-iso")) %>%
    tm_map(removeWords, filterWords) %>% 
    tm_map(FUN=removePunctuation) %>%
    tm_map(FUN=stripWhitespace)
if (prod(is.na(additionalRemove))!=1){
    docs <- docs %>% tm_map(removeWords, additionalRemove)
}
docs <- TermDocumentMatrix(docs)
mat <- as.matrix(docs)
matSorted <- sort(rowSums(mat), decreasing=TRUE)
    
## Perform the same filtering for whole data
allDocs <- VCorpus(VectorSource(tb$Gene_summary))
allDocs <- allDocs %>%
    tm_map(FUN=content_transformer(tolower)) %>% 
    tm_map(FUN=removeNumbers) %>%
    tm_map(removeWords, stopwords::stopwords("english", "stopwords-iso")) %>%
    tm_map(removeWords, filterWords) %>% 
    tm_map(FUN=removePunctuation) %>%
    tm_map(FUN=stripWhitespace)
if (prod(is.na(additionalRemove))!=1){
    allDocs <- allDocs %>% tm_map(removeWords, additionalRemove)
}
matAll <- as.matrix(TermDocumentMatrix(allDocs))
matAllSorted <- sort(rowSums(matAll), decreasing=TRUE)
        

## ORA
returnP <- function(name){
    query <- as.numeric(matSorted[name])
    noquery <- sum(matSorted) - query
    queryAll <- as.numeric(matAllSorted[name])
    allwords <- sum(matAllSorted) - queryAll

    ## p-value
    return(sum(dhyper(query:sum(matSorted), queryAll, allwords, sum(matSorted))))
}

sigs <- names(matSorted)[p.adjust(sapply(names(matSorted), returnP), "bonferroni")<0.05]
sigs
