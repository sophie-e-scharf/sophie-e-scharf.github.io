---
title: "Signal Detection Model with STAN - Part 1"
author: "David Izydorczyk"
date: 2020-04-28
output: 
  html_document:
    keep_md: true
    toc: true
    df_print: paged
    number_sections: true
    theme: paper
    highlight: tango
bibliography: "../Literatur/my_bib.bib"
csl: "../apa6.csl"
layout: post
categories: R
# status: process
# published: true
status: publish
published: true
---
 
 
- Describe Model + Parameters more
- Describe results more
 
# About this blog post
 
This blog post is supposed to be the first blog post in a series of blog posts comparing the SDT model with a 2THM model. I wrote this blog post for three reasons. First, I wanted to learn more about the modeling of memory data. Second, I wanted to learn more about STAN, since I mostly use JAGS in my own research. Finally, I also wanted to get more practice in writing, since I take forever when writing my own articles (which I am supposed to do right now, when writing this blogpost). So the reasons for this (and following blog posts) are rather selfish. However, if anyone ever finds this blog post and finds it helpful, that would be even better !
 
This blog post is about modeling data from a memory recogntion task with a non-hierarchical SDT model. In following blog posts, I will extend this model to account for differences between individuals as well as differences between stimulus sets. Then I will model the same data with a 2HTM and finally compare both models with each other. 
 
# Setup 
 
At the beginning I the load packages I need for this analysis and set the general settings for the code chunks of the markdown document.
 

{% highlight r %}
library(tidyverse) # contains ggplot, dplyr etc
library(rstan)
library(bayesplot)
library(rstanarm)
library(patchwork)
library(truncnorm)
library(kableExtra)
library(knitr)
{% endhighlight %}
 

{% highlight r %}
knitr::opts_chunk$set(echo  = TRUE,cache = TRUE)
bayesplot_theme_set(theme_bw())
{% endhighlight %}
 
# Signal detection model
 
In this part, I want to model the data with a signal detection model (SDT) as described in Chapter 11 of the fantastic book from @lee2014. Signal detection theory (SDT) may be applied to any area of psychology in which two different types
of stimuli must be discriminated. It was first applied in studies of perception, where subjects had to discriminated between signals (stimuli) and noise (no stimuli). However, SDT is also often used to decribe  recognition memory, where subjects have to discriminate between old (signal) and new items (noise). 
 
The idea behind SDT is that signal and noise trials can be represented as values along a uni-dimensional "strength" dimension, where signal trials are assumed to have a greater strength than noise trials. According to SDT, people produce "old" or "new" decisions by comparing the strength of the current trial to a fixed threshold. If the strength exceeds this threshold, the response is "old", otherwise the response is "new". The distributions of strength of signal and noise trials are assumed to be Normal distributions with different means (the mean of the noise distributions is equal to 0), but the same variance (which is fixed to 1), which can be expressed as:
 
$$
noise \sim N(0,1) \\
signal \sim N(d,1)
$$
where *d* is the *discriminability* or *sensitivity* parameter, which goes from $-\infty$ to $+ \infty $. This parameter *d* corresponds to the distance between means of the noise and the signal distribution in standard deviation units. A value of 0 indicates an inability to distinguish signals from noise, whereas larger values indicate a correspondingly greater ability to distinguish signals from noise. Negative values are also possible, but are harder to interpret. They are most often thought of as a sampling error or response confusion (i.e., responding "old" when intending to respond "new", and vice vers, for instance by confusing the corresponding buttons). 
 
Another parameter is the *c*, the *bias* parameter. Positive values of *c* correspond to a bias towards saying *new* and negative values correspond to a bias towards saying *old*. 
 
These two parameters can then be directly translated into hit (HR) and false-alarm rate (FR) via:
 
$$
HR = \phi(0.5 \times d - c) \\
FR = \phi(- 0.5 \times d - c) 
$$
 
where, $\phi$ is the cumulative density function of the standard normal distribution. HR and FR  map naturally to the data pattern we observe in a recognition memory experiment: 
 

|Response |Signal Trial |Noise Trial |
|:--------|:------------|:-----------|
|Old      |Old          |New         |
|New      |Old          |New         |
 
which corresponds to:
 

|Response |Signal Trial |Noise Trial       |
|:--------|:------------|:-----------------|
|Old      |Hit          |False alarm       |
|New      |Miss         |correct rejection |
 
 
 
This allows us to take the observed number of hits and false alarms in our experiment, and translate them into psychological meaningfull parameters *d* and *c*
 
# The Experiment 
 
The experiment was conducted as part of a pre-test to select stimuli for a subsequent multiple-cue judgment experiment. The experiment consisted of five blocks, with two phases each. In the learning phase, participants saw 12 different pictures of either plants, houses, aliens, bugs, or plates, (a different stimulus set in each block) with features differing on five cues. The same 12 pictures were presented 4 times for 5 seconds to each participant. After the learning phase, participants were presented again with all 12 old pictures as well as 20 new pictures. They were asked to indicate if the picture was an old picture or a new one.  These two phases were repeated for each of the five stimulus sets. 
 
# The Data set 
 
The data I will use contains the following seven variables:
 
- `ID`: participants ID
- `block`: block number 
- `trial`: trial number
- `stimulus`: type of stimulus (house, alien, plant, plate, or bug)
- `response`: response given in this trial (old or new)
- `correctResponse`: correct response in this trial
- `correct`: is `true` if response is equal to `correctResponse`, if not is it `false`
 

{% highlight r %}
data <- read.csv2("../Data/Stimulus_Test_tidy.csv") %>% 
            select(ID,block,trial,stimulus,response,correctResponse,correct)
{% endhighlight %}



{% highlight text %}
## Warning in file(file, "rt"): cannot open file '../Data/
## Stimulus_Test_tidy.csv': No such file or directory
{% endhighlight %}



{% highlight text %}
## Error in file(file, "rt"): cannot open the connection
{% endhighlight %}



{% highlight r %}
data[1:6,]
{% endhighlight %}



{% highlight text %}
## Error in data[1:6, ]: object of type 'closure' is not subsettable
{% endhighlight %}
 
Based on these variables, I calculated the number of hits, false alarms, false negatives, and misses, as well as their corresponding rates for each person, in each block. In addition, I calculated $d'$  as $d' = z(h) - z(fa)$ [@snodgrass1988; @stanislaw1999], where $h$ is the hit rate and $fa$ is the false-alarm rate.
 

{% highlight r %}
hits <- data %>%  
          mutate(IDn            = as.numeric(ID), 
                 hit            = ifelse(response == "old" & correctResponse == "old",1,0),
                 false_positive = ifelse(response == "old" & correctResponse == "new",1,0)) %>% 
          group_by(ID,IDn,stimulus,block) %>% 
          summarize(n_old    = sum(correctResponse == "old"),
                    n_new    = sum(correctResponse == "new"),
                    h        = sum(hit),
                    fa       = sum(false_positive),
                    hit_rate = sum(hit)/n_old, 
                    fa_rate  = sum(false_positive)/n_new) %>% 
          mutate(hit_rate = case_when(
                               hit_rate == 1 ~ .9999999,
                               hit_rate == 0 ~ .0000001,
                               TRUE ~ hit_rate
                            ),
                 fa_rate = case_when(
                               fa_rate == 1 ~ .9999999,
                               fa_rate == 0 ~ .0000001,
                               TRUE ~ fa_rate
                            ),
                 z_h_rate  = qnorm(hit_rate),
                 z_fa_rate = qnorm(fa_rate),
                 dprime    = z_h_rate-z_fa_rate)  
{% endhighlight %}



{% highlight text %}
## Error in UseMethod("mutate_"): no applicable method for 'mutate_' applied to an object of class "function"
{% endhighlight %}
 
Which then gives us the following data structure:
 

{% highlight r %}
hits[1:5,]
{% endhighlight %}



{% highlight text %}
## Error in eval(expr, envir, enclos): object 'hits' not found
{% endhighlight %}
 
 
I will first build the non-hierarchical version of the SDT ^[Described on pages 158-159 in @lee2014] analyzing the data from only one stimulus set. 
 
## SDT - Non Hierarchical - One Stimulus Set 
 
As a beginning, I will select the data from only specific stimulus set. 
 

{% highlight r %}
temp <- hits %>% filter(stimulus == "plants")
{% endhighlight %}



{% highlight text %}
## Error in eval(lhs, parent, parent): object 'hits' not found
{% endhighlight %}
 
 
### The Model 
 
Next, I will define the SDT model in STAN.  
 
#### The Data 
 
In the `data` block, I define the data we are using in the model, which are:
 
- `p` the number of participants
- `s` a vector containing the number of signal trials (i.e., trials with a picture already shown in the learning phase)
- `n` a vector containing the number of noise trials (i.e., the number of new pictures)
- `h` a vector containing the number of hits of each person
- `fa` a vector containing the number if false alarms of each person
 
 

{% highlight r %}
stan_data <- list(
  s  = temp$n_old,
  n  = temp$n_new,
  h  = temp$h,
  fa = temp$fa,
  p  = length(unique(temp$ID))
)
{% endhighlight %}



{% highlight text %}
## Error in eval(expr, envir, enclos): object 'temp' not found
{% endhighlight %}



{% highlight r %}
stan_data
{% endhighlight %}



{% highlight text %}
## Error in eval(expr, envir, enclos): object 'stan_data' not found
{% endhighlight %}
 
#### The Parameters 
 
In the `parameters` block, I define the two main parameters we have in our model and which we are interested in, *d* and *c*. These two parameters are then transformed in to the *hit rate* and the *false alarm rate* in the `transformed parameters` block. In the `model` block I then define the binomial likelihood function connecting our parameters and our data. In the `generated quantities` block I also include variables for later posterior predictive analysis, capturing the predictinos of `h` and `fa` based on the current parameter values of each step of the MCMC-chains.
 
#### The Priors 
 
The priors for both *d* and *c* are normal distributions with $\mu = 0$ and $\sigma = 1$, which corresponds to a uniform distribution after transforming into hit and false-alarm rates.
 
 
![plot of chunk unnamed-chunk-4](/assets/img/SDT_nonhierarchical.Rmd/unnamed-chunk-4-1.png)
 
#### The STAN-Model 
 
The model code then looks like this:
 

{% highlight stan %}
data {
  int<lower=0> p;     // number of persons
  int<lower=0> s[p];  // number of signal trials == old trials of person p
  int<lower=0> n[p];  // number of noise  trials == new trials of person p
  int<lower=0> h[p];  // number of hits of person p
  int<lower=0> fa[p]; // number of false alarms of person p
}
parameters {
  vector[p] d;
  vector[p] c;
}
transformed parameters {
  vector<lower=0,upper=1> [p] hit_rate;
  vector<lower=0,upper=1> [p] fa_rate;
 
  for (i in 1:p){
    hit_rate[i] = Phi_approx(d[i]/2-c[i]);
    fa_rate[i]  = Phi_approx(-d[i]/2-c[i]);
  }
}
model {
// Priors
  d ~ normal(0, 1);
  c ~ normal(0, 1);
 
// Likelihood
  h  ~ binomial(s, hit_rate);
  fa ~ binomial(n, fa_rate);
} generated quantities {
 
  vector<lower=0> [p] h_pred;
  vector<lower=0> [p] fa_pred;
 
  for (i in 1:p){
    h_pred[i]  = binomial_rng(s[i], hit_rate[i]);
    fa_pred[i] = binomial_rng(n[i], fa_rate[i]);
  }
}
{% endhighlight %}
 
### Run the Model 
 
Next, we sample from the model using 10.000 iterations, with  a rather small warm-up of 2000 iterations and a thin = 4. 
 

{% highlight r %}
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
 
fit1 <-sampling(
        object     = StanModel_oneStim_nonHier, # Stan program
        data       = stan_data, # named list of data
        chains     = 2,         # number of Markov chains
        warmup     = 2000,      # number of warmup iterations per chain
        iter       = 10000,     # total number of iterations per chain
        cores      = 2,         # number of cores (could use one per chain)
        refresh    = 0,         # no progress shown
        thin       = 4
      )
{% endhighlight %}



{% highlight text %}
## Error in .local(object, ...): object 'stan_data' not found
{% endhighlight %}



{% highlight r %}
# save results 
posterior_SDT_oS_nH <- rstan::extract(fit1, permuted = FALSE)
{% endhighlight %}



{% highlight text %}
## Error in rstan::extract(fit1, permuted = FALSE): object 'fit1' not found
{% endhighlight %}



{% highlight r %}
summary_oS_nH       <- summary(fit1) %>% as.data.frame()
{% endhighlight %}



{% highlight text %}
## Error in summary(fit1): object 'fit1' not found
{% endhighlight %}
 
### Inspect MCMC, Rhat, ESS
 
First, we can look at some MCMC-Traces for some of the parameters and persons.
 

{% highlight r %}
 mcmc_trace(posterior_SDT_oS_nH,
            pars = vars("d[1]":"d[3]",
                        "c[1]":"c[3]",
                        "hit_rate[1]":"hit_rate[3]",
                        "fa_rate[1]" :"fa_rate[3]"),
            facet_args = list(nrow = 4, labeller = label_parsed))
{% endhighlight %}



{% highlight text %}
## Error in is.data.frame(x): object 'posterior_SDT_oS_nH' not found
{% endhighlight %}
 
So far so good, this looks exactly as you would like it, some nice hairy caterpillars.
 
 
We can also look at the distributions of *Rhats*  and the effective sample size:
 

{% highlight r %}
p1 <- ggplot(summary_oS_nH,aes(x =summary.n_eff))+
        geom_histogram(bins = 40, color = "black", fill = "skyblue2")+
        theme_bw() +
        labs(x = "Effective Sample Size",
             y = "Count")
{% endhighlight %}



{% highlight text %}
## Error in ggplot(summary_oS_nH, aes(x = summary.n_eff)): object 'summary_oS_nH' not found
{% endhighlight %}



{% highlight r %}
p2 <- ggplot(summary_oS_nH,aes(x =summary.Rhat))+
        geom_histogram(bins = 40, color = "black", fill = "tomato2")+
        theme_bw() +
        labs(x = "Rhat",
             y = "Count")
{% endhighlight %}



{% highlight text %}
## Error in ggplot(summary_oS_nH, aes(x = summary.Rhat)): object 'summary_oS_nH' not found
{% endhighlight %}



{% highlight r %}
p1+p2
{% endhighlight %}

![plot of chunk unnamed-chunk-8](/assets/img/SDT_nonhierarchical.Rmd/unnamed-chunk-8-1.png)
 
This looks also good. The effective sample sizes are always > 400, which is often recommended (citation), and are near the total sample size of  4000.
 
### Posterior Summaries & Densities
 
Next, we can look at the summary statistics and plots of the posterior densitie distributions. So lets make a convient tidy data.frame for plotting and for the posterior summary statistics. 
 

{% highlight r %}
desc_oS_nH_tidy <- summary_oS_nH %>% 
                      mutate(rn     = rownames(.),
                              param = str_split(rn,"\\[", simplify = TRUE) %>% unlist() %>%  .[,1],
                              ID    = str_extract(rn,"\\d+") %>% as.numeric()) %>% 
                      select(param,ID,
                             mean    = summary.mean,
                             se_mean = summary.se_mean,
                             HDI025  = summary.2.5.,
                             median  = summary.50.,
                             HDI975  = summary.97.5.) 
{% endhighlight %}



{% highlight text %}
## Error in eval(lhs, parent, parent): object 'summary_oS_nH' not found
{% endhighlight %}



{% highlight r %}
mcmc_SDT_oS_nH <-  rstan::extract(fit1,pars="lp__",include=FALSE) %>%
                      bind_rows() %>% 
                      mutate(ID     = rep(1:40,each = 4000),
                             dprime = rep(temp$dprime, each = 4000)) %>% 
                      group_by(ID) %>% 
                      mutate(mPP_d  = median(d),
                             mPP_c  = median(c))
{% endhighlight %}



{% highlight text %}
## Error in rstan::extract(fit1, pars = "lp__", include = FALSE): object 'fit1' not found
{% endhighlight %}
 
 
#### d 
 
Since *d* is the parameter we are most interested in, lets start with this one. Below are the descriptive statistics, a forest plot showing the median, 50 \%, and 96 \% credible intervals of the posterior distributions, as well as univariate marginal posterior distributions for some participants, showing the median posterior estimate (red dashed line) as well as the analyticaly calculated *d'* (black dashed line).
 

{% highlight r %}
desc_oS_nH_tidy %>% filter(param %in% c("d"))
{% endhighlight %}



{% highlight text %}
## Error in eval(lhs, parent, parent): object 'desc_oS_nH_tidy' not found
{% endhighlight %}
 
 

{% highlight r %}
mcmc_intervals(posterior_SDT_oS_nH,
               pars = vars(starts_with("d")),
               point_est = "median",
               prob = 0.5,
               prob_outer = 0.96)
{% endhighlight %}



{% highlight text %}
## Error in is.data.frame(x): object 'posterior_SDT_oS_nH' not found
{% endhighlight %}
 
 

{% highlight r %}
mcmc_SDT_oS_nH %>%
  filter(ID %in% sample(1:40,size=4)) %>% 
  ggplot(., aes(x = d))+
    geom_density(color = "black",fill = "skyblue2",alpha=0.5) +
    geom_vline(aes(xintercept = mPP_d), color = "red", lty = "dashed",lwd = 1) +
    geom_vline(aes(xintercept = dprime), color = "black", lty = "dashed") +
    facet_grid(.~ID,scales="free") + 
    theme_bw()
{% endhighlight %}



{% highlight text %}
## Error in eval(lhs, parent, parent): object 'mcmc_SDT_oS_nH' not found
{% endhighlight %}
 
We can see that the *d'* values are rather small and very close to 0 for most people. This indicates that participants had a hard time differentiating old from new stimuli. I had hoped for larger values, as it is necessary that people are able to discriminate between stimuli and their feature rather well for the multiple-cue judgment experiment I want to conduct with these stimuli.
 
Also, from the plots it is also evident that the median of the posterior distribution is very close to the analytically calculated *d'* value.  
 
 
#### c
 
We can also look at the *c* values. Since we have more noise trials (new stimuli) than signal trials (old stimuli) in our testing phase (20 vs. 12), and participant were told this information, I expected to find slightly positive values of *c*. However, as apparent from the summary statistics and the plots, participants had more negative values of *c*, indicating a bias for the "old"-response.
 

{% highlight r %}
desc_oS_nH_tidy %>% filter(param %in% c("c"))
{% endhighlight %}



{% highlight text %}
## Error in eval(lhs, parent, parent): object 'desc_oS_nH_tidy' not found
{% endhighlight %}
 

{% highlight r %}
mcmc_intervals(posterior_SDT_oS_nH,
               pars = vars(starts_with("c")),
               point_est = "median",
               prob = 0.5,
               prob_outer = 0.96)
{% endhighlight %}



{% highlight text %}
## Error in is.data.frame(x): object 'posterior_SDT_oS_nH' not found
{% endhighlight %}
 

{% highlight r %}
mcmc_SDT_oS_nH %>%
  filter(ID %in% sample(1:40,size=4)) %>% 
  ggplot(., aes(x = c))+
    geom_density(color = "black",fill = "skyblue2",alpha=0.5) +
    geom_vline(aes(xintercept = mPP_c), color = "red", lty = "dashed",lwd = 1) +
    facet_wrap(.~ID,scales="free") + 
    theme_bw()
{% endhighlight %}



{% highlight text %}
## Error in eval(lhs, parent, parent): object 'mcmc_SDT_oS_nH' not found
{% endhighlight %}
 
### The posterior predictive values
 
We can also look at the posterior predictive values of `h` and `fa`. 
 

{% highlight r %}
temp     <- hits  %>%  filter(stimulus == "plants") %>% ungroup() %>% select(.,IDn,h,fa)
{% endhighlight %}



{% highlight text %}
## Error in eval(lhs, parent, parent): object 'hits' not found
{% endhighlight %}



{% highlight r %}
postpred <- mcmc_SDT_oS_nH %>% 
              select(., ID, h_pred, fa_pred) %>% 
              left_join(.,temp,by = c("ID" = "IDn")) 
{% endhighlight %}



{% highlight text %}
## Error in eval(lhs, parent, parent): object 'mcmc_SDT_oS_nH' not found
{% endhighlight %}



{% highlight r %}
head(postpred)
{% endhighlight %}



{% highlight text %}
## Error in head(postpred): object 'postpred' not found
{% endhighlight %}
 
 

{% highlight r %}
pp1 <- ggplot(filter(postpred,ID %in% 1:4),aes(x = h_pred)) +
        geom_histogram(bins=10,color = "black",fill = "skyblue2",alpha=0.5) +
        geom_vline(aes(xintercept =  h),color = "black",lwd = 1.5, lty = "dashed")+
        facet_wrap(.~ID,scales="free") + 
        theme_bw() 
{% endhighlight %}



{% highlight text %}
## Error in filter(postpred, ID %in% 1:4): object 'postpred' not found
{% endhighlight %}



{% highlight r %}
pp2 <- ggplot(filter(postpred,ID %in% 1:4),aes(x = fa_pred)) +
        geom_histogram(bins=10,color = "black",fill = "tomato2",alpha=0.5) +
        geom_vline(aes(xintercept =  fa),color = "black",lwd = 1.5, lty = "dashed")+
        facet_wrap(.~ID,scales="free") + 
        theme_bw() 
{% endhighlight %}



{% highlight text %}
## Error in filter(postpred, ID %in% 1:4): object 'postpred' not found
{% endhighlight %}



{% highlight r %}
pp1 + pp2
{% endhighlight %}



{% highlight text %}
## Error in eval(expr, envir, enclos): object 'pp1' not found
{% endhighlight %}
 
And we can also calculate how often the 95% credible interval of the predicted values contains the true values:
 

{% highlight r %}
temp1 <- hits  %>%  filter(stimulus == "plants") %>% ungroup() %>% select(.,ID=IDn,h,fa)
{% endhighlight %}



{% highlight text %}
## Error in eval(lhs, parent, parent): object 'hits' not found
{% endhighlight %}



{% highlight r %}
temp2 <- desc_oS_nH_tidy %>% 
          select(param,ID,HDI025,HDI975)  %>% 
          filter(param == "h_pred" | param == "fa_pred") %>% 
          pivot_wider(.,id_cols    = ID,
                        names_from = param,
                                          values_from = c(HDI025,HDI975))
{% endhighlight %}



{% highlight text %}
## Error in eval(lhs, parent, parent): object 'desc_oS_nH_tidy' not found
{% endhighlight %}



{% highlight r %}
left_join(temp1,temp2,by = "ID")  %>% 
  summarize(pp_h  = mean(h > HDI025_h_pred & h < HDI975_h_pred),
            pp_fa = mean(fa > HDI025_fa_pred & fa < HDI975_fa_pred))
{% endhighlight %}



{% highlight text %}
## Error in left_join(temp1, temp2, by = "ID"): object 'temp1' not found
{% endhighlight %}
 
As evident from the plots and the 95% credible interval checks of the posterior predictives, our model is able to recover our data well.
 
