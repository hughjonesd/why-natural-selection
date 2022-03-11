
# NHB reviewer comments to deal with

## Rev 1

3) The integration of an economic theory explaining the observed patterns is in general interesting and desirable, but a) would have been more interesting when established in advance, b) has ambiguous predictions and appears to be utilised post-hoc and then not clearly tested, c) has nothing to do with genetic research, but is a behavioural approach and it thus d) would be very surprised if it is not yet established in some way in the non-genetic literature and tested in such a context as the UK birth cohorts studies. Yet, I have not seen any related demographic or economic literature been criticised.

- Just wrong about the related dem/econ lit. Can point out
  the theory is old. Clarify the testing.
  
4) Regarding potential implications on inequality, e. g. in health, maybe a quantification would have been interesting but the way it is presented, it is too superficial to be a clear advance in the field.

- Could we quantify? E.g. if genetics correlation with the phenotype (h2) is x,   then what does an increase in variance of genetics imply about increase in
  variance of the phenotype?

Weighted results differ up to a factor of 1.8, and I didn’t understand how the authors deal with this. 

- Clarify where we use weights and robustness. CHECK

The inclusion of two generations would be a particularly exciting feature, but the fact that the first generation is selected for having a child does significantly reduce this benefit. The authors also don’t appear to deal with this issue by e. g. testing for linearity in effects by parity or applying some sort of selection model, so this remains an important issue.

- We could test for linearity in effects. I think our approach is still preferable.

he sample of 400k looks quite large even for UKB and I believe that individuals have been selected to be of European ancestry, but not all are born in the UK. I would suggest only to use such individuals to reduce selection and avoid further heterogeneity.

- We could do. Unlikely to reduce selection. I'd say, meh.

I do not understand this statement:
“The “selection effect,” 𝛽, reflects the strength of natural selection within the sample. In fact, since polygenic scores are normalized, 𝛽 is the expected polygenic score among children of the sample”. Is it the expected value of the polygenic score amongst the children, is it the average? Is it the expected change of the PGS?

- Clarify language in footnote. CHECK. The footnote is fine, guy obviously too lazy to read.

explain what balancing and diversifying selection is and how you test it exactly?

- CHECK language. 
- DONE: altered language in appendix.

I do not understand how you exactly understand reasons for ascertainment bias and learn about mechanisms for natural selection at the same time by splitting the sample by SES

- CHECK language.
- DONE. Rewrote to focus on motivating theory.

For the split data in Figure 3: do you use weights?

- Maybe clarify weighting throughout. CHECK.

## Rev 2

The authors presentation of their results needs improvement. There are so many results, plus a model, and as a result the paper reads like a disconnected set of findings and lacks an overall coherent narrative. In my opinion, the paper’s most important finding is that human capital mediates natural selection. I suspect the authors agree with me on this point, since the title of their manuscript says just that. Still, when reading the paper, one does not get to the topic of human capital and natural selection until p. 11 (and the paper has only 16 pages)! I would suggest something along these lines:

> After the introduction, dive directly into the topic of human capital and natural selection, and devote more space to that topic (include Figure 1 in there). For instance, Figure 6 should be shown early on, and Appendix Figure 17 should go in the main text. I also have other requests or suggestions in point 2 below.

- CHECK and consider. Risk is then we pretend theory came first. 
- DONE to some extent in that we're putting theory earlier, and clarifying the motivation for the subgroups.

There’s a disconnect between the lengthy discussion of the stratified analysis (Figures 3-5 and other Appendix figures) – which ends by concluding that AFLB mediates the association between the PGIs and number of children – and the title of the paper and the model, which are about human capital. That long discussion is interesting but also distracts from what I think should be the main point of the paper. Plus, the stratified analyses are limited and only so much can be learned from them. As the authors acknowledge at lines 73-74: “this exercise does not apportion the total correlation with fertility into subgroups. That is because the polygenic score may also affect the probability of being in each subgroup.” Thus, I recommend that that discussion be substantially shortened and presented towards the end of the paper, with most of the figures and details relegated to the Appendix.

- Ditto. This is just because they don't like it. "Please hide these nasty results away."

given that the authors’ claims are about mediation, it’d be nice to for them to also do a full-fledged mediation analysis and to estimate the share of each PGS-Number of children association that’s mediated by the different measures of human capital together (not only one by one). For an example of such a framework, the authors could refer to Section 6 of the Supplementary Information of the “EA2” paper (Okbay et al., Nature Genetics 2016). 

- Probably "full-fledged mediation analysis"" isn't, but perhaps Oana's framework could be reused to say under what conditions etc.
- Indeed Okbay have a (thoughtful & modest) discussion of the theory behind their mediation analysis
  which is indeed quite standard. 
- Could we use within-siblings to estimate PGS on educ/earnings? 
  And then within-siblings to estimate PGS on fertility?
  And then get estimates for earnings/educ on fertility from somewhere else "causal"?
  Or just regress earnings/educ on fertility in the standard vulgar way?
  Seems like the same approach that we use in the paper with Oana
  The Okbay appendix explains how to calculate stderrs via delta method

To adjust for ascertainment bias, the authors adjust their estimates using three weighting schemes. Why not combine all the variables from the three schemes and add other observables that could impact ascertainment (including also sex, proximity to a recruitment center, household income)? 

- Maybe get weights from the other guys. [Sent email]

Also, instead of presenting a small subset of the results with various weighting schemes in Figure 2, it’d be preferable for the authors to pick what they think is the best weighting scheme and then present all results with that scheme. Then, they can have a paragraph talking about the robustness of their results to various weighting schemes and show some of these checks in the Appendix.

- That might make sense.

How were the 33 PGSs selected? 

- We don't have a real answer but hey... 

Also, they should mention another source of bias vs. the analysis with the children: errors-in-variables in their measure of the parents’ scores (which is taken as the average of the children’s scores, as opposed to being measured directly in the parents).

- Think about. Not sure this is big deal because both parents indeed should
average children's scores closely.
- No, there could be error, simply by randomness of meiosis. DONE.

12. In the Results section, at around line 44, the authors should specify that all 409,629 individuals are of European descent (they say that in Materials and methods, but this should also be mentioned in the main text).

- Sure. DONE

# Rev 3

 That selection appears strongest in the youngest population has other problems: heritability estimates often increase with time. The point being: the selection effects are strongest specifically in the population that is most dependent on their environment for resources, and consequently, where we should trust the polygenic scores the least (and in particular for behavioral and cognitive traits). 
 
 - This doesn't make sense. If you trust the PGS at age 30, then they
 are also interesting at age 20, b/c 20 year olds become 30 year olds.
 
Everything else is just inane bollocks, except we ought to lay out the
2 alt theories we test, and add the theory of welfare to point out it
isn't really compatible with no-change-over-time.

DONE



