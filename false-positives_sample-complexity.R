

# ===============================
# basically just copied this code from the false negative analysis
#
# one note: for the 'pollen.grain.proportion' analysis the data are aggregated by mix
# see the section 'format kraken data for false positive analysis'
# I used the *minimum* pollen grain proportion for any taxon in a mix as the value here
# ... so the interpretation will be slightly different than for the false negative analysis

# ===============================

# set up the three factors by which we are running the models
taxon = c("species", "genus", "family")
datasubset = c("sub", "all") # whether we are using the designated subset of data designed for the question, or all data
question = c("spp.rich", "relatedness", "pollen.grain.proportion")

# calculate total number of models
total = length(taxon)*length(datasubset)*length(question)

# first set up a table for the results with a number of entries equal to the 'total' variable above (18):
results.table.falsepos = data.frame(question = rep(NA,total), taxon  = rep(NA,total), data.subset  = rep(NA,total), model.name = rep(NA,total), p.val  = rep(1.000001,total), n  = rep(9999,total), warning.msg = rep(NA,total))

# keep track of which row of the table to record in:
tracker = 1

# EXAMPLE FORMULA:
# Krak.Q1.species.all = glmer(qual.species.rbcL ~ spp.rich + (1|mix.ID/sample/rep.ID) + (1|species), family = binomial, data = agg.krak.species, control = glmerControl(optimizer="bobyqa"))

# 'for' loop:

for(q in 1:3) { # 'question': response variables for Q1 / Q2 / Q3
  for(k in 1:3){ # 'taxon': species, genus, family
    for(l in 1:2) { # 'datasubset': sub or all
      # first, name the analysis:
      namer = paste("Krak.falsepos.Q", q, ".", taxon[k], ".", datasubset[l], sep = "")
      # second, set which taxonomic data to use:
      data.to.use = paste("agg.krak.", taxon[k], sep = "")
      # third, set up the data subset
      subster = paste("data.sub = filter(", data.to.use, ", question.1 == ", q, " | question.2 == ", q, " | question.3 == ",q, ")", sep = "")
      eval(parse(text = subster)) # probably not the most efficient thing ever... 
      # fourth, set whether or not data subset is used (vs. all data)
      if(datasubset[l]=="sub") {data.to.use = "data.sub"} # i.e., doesn't change if all data are to be used
      # fifth, set up mixed-effects model: 
      mixed = paste(namer, " = suppressWarnings(glmer(qual.", taxon[k],
                    " ~ ",  question[q], " + (1|mix.ID/krak.sample/krak.rep.ID) + (1|", taxon[k], "), family = binomial, 
                    data = ", data.to.use, ", control = glmerControl(optimizer=\"bobyqa\")))", sep = "")
      # sixth, evaluate the mixed-effects model
      eval(parse(text = mixed))
      # # eighth, print summary of model [SKIP FOR NOW]
      # summarizer = paste("print(summary(", namer, "))", sep = "")
      # eval(parse(text = summarizer)) # print summary of the mixed-effects model
      
      ## extract p-value
      # (this would probably be easier using the 'broom.mixed' package?)
      # example: coef(summary(Q3.genus.ITS.all))[2,4]
      pvaller = paste("pval <- coef(summary(", namer, "))[2,4]", sep = "")
      eval(parse(text = pvaller))
      
      # extract convergence failures
      converger = paste(namer, "@optinfo$conv$lme4$code", sep = "")
      converg = eval(parse(text = converger))
      converg.return = ifelse(length(converg)==1, "ERROR!!", "")
      
      # record results in table
      results.table.falsepos[tracker,1] = question[q]
      results.table.falsepos[tracker,2] = taxon[k]
      results.table.falsepos[tracker,3] = datasubset[l]
      results.table.falsepos[tracker,4] = namer
      results.table.falsepos[tracker,5] = pval
      results.table.falsepos[tracker,6] = nrow(eval(parse(text = data.to.use)))
      results.table.falsepos[tracker,7] = converg.return
      
      # advance tracker
      tracker = tracker + 1
    }
  }
}

# display results table
# note that the 'kable' function is part of the 'knitr' package and `kable_styling` is from the `kableExtra` package 
kable(results.table.falsepos) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F)

