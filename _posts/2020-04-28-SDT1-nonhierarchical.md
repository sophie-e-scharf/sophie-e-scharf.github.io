---
layout: single
title: "Signal Detection Model with STAN - Part 1"
date: 2020-04-28
tags:
  - cognitive modeling
  - STAN
  - R
  
# status: process
# published: true
status: publish
published: true
permalink: /posts/2020/04/SDT-1/
---
 

 
 
 
## About this blog post
 
This blog post is supposed to be the first in a series of blog post comparing the *signal detection theory * model (SDT) with a *two-high-threshold* model (2HTM) of recognition. I wrote this blog post for three reasons. First, I wanted to learn more about the modeling of recognition and memory data. Second, I wanted to learn more about STAN, since I mostly use JAGS in my own research. Finally, I also wanted to practice writing, since I take forever when writing my own articles. So the reasons for this blog are rather selfish. However, if anyone ever finds this blog post and finds it helpful, that would be even better!
 
 
In this blog post you are going to read, I will use a *non-hierarchical SDT* model to investigate the data from a recogntion experiment. In following blog posts, I will extend this model to account for differences between individuals as well as differences between stimulus sets. Then I will model the same data with a 2HTM and finally compare both models with each other. 
 
## Setup 
 
At the beginning, I load the packages I need for this analysis. However, before we start with the actual analysis and modeling, I first want to give a short introduction into signal detection theory.
 

{% highlight r %}
library(tidyverse) # contains ggplot, dplyr etc
library(rstan)
library(bayesplot)
library(rstanarm)
library(patchwork)
library(truncnorm)
library(kableExtra)
library(knitr)
 
# Set bayesplot theme
bayesplot_theme_set(theme_bw())
{% endhighlight %}
 
 
## Signal detection model
 
Signal detection theory (SDT) may be applied to any area of psychology in which two different types of stimuli must be discriminated. It was first applied in studies of perception, where subjects had to discriminated between signals (stimuli) and noise (no stimuli). However, SDT is also often used to describe memory recognition task, where subjects have to discriminate between old (signal) and new items (noise). 
 
The idea behind SDT is that signal and noise trials can be represented as values along a uni-dimensional continuous strength dimension, where signal trials are assumed to have a greater strength than noise trials. According to SDT, people produce "old" or "new" decisions by comparing the strength of the current trial to a fixed threshold. If the strength exceeds this threshold, the response is "old", otherwise the response is "new". The strength of signal and noise trials is assumed to be normally distributed with different means (but the mean of the noise distributions is assumed to be equal to 0), but the same variance (which is fixed to 1), which can be expressed as:
 
$$
\begin{aligned}
noise &\sim N(0,1)  \\[.5em]
signal &\sim N(d,1)
\end{aligned}
$$
 
where **d** is the **discriminability** or **sensitivity** parameter, which goes from $-\infty$ to $+ \infty $. This parameter **d** corresponds to the distance between means of the noise and the signal distribution in standard deviation units. A value of 0 indicates an inability to distinguish signals from noise, whereas larger values indicate a correspondingly greater ability to distinguish signals from noise. Negative values are also possible, but are harder to interpret. They are most often thought of as a sampling error or response confusion (i.e., responding "old" when intending to respond "new", and vice versa, for instance by confusing the corresponding buttons). 
 
Another parameter is **c**, the **bias** parameter. Positive values of **c** correspond to a bias towards saying *new* and negative values correspond to a bias towards saying *old*. 
 
These two parameters can then be directly translated into hit (HR) and false-alarm rates (FR) via:
 
$$
\begin{aligned}
HR &= \phi(0.5 \times d - c) \\[.5em]
FR &= \phi(- 0.5 \times d - c) 
\end{aligned}
$$
 
where, $\phi$ is the cumulative density function of the standard normal distribution. HR and FR  map naturally to the data pattern we observe in a recognition memory experiment, were participants see either an old (signal trial) or a new item (noise trial), and then have to respond by pressing the corresponding button:
 

|Response |Signal Trial |Noise Trial |
|:--------|:------------|:-----------|
|Old      |Old          |New         |
|New      |Old          |New         |
 
which corresponds to:
 

|Response |Signal Trial |Noise Trial       |
|:--------|:------------|:-----------------|
|Old      |Hit          |False alarm       |
|New      |Miss         |correct rejection |
 
This allows us to take the observed number of hits and false alarms in our experiment, and translate them into psychological meaningful parameters *d* and *c*
 
## The Experiment 
 
The experiment was conducted as part of a pre-test to select stimuli for a subsequent multiple-cue judgment experiment. The experiment consisted of five blocks, with two phases each. In the learning phase, participants saw 12 different pictures of either plants, houses, aliens, bugs, or plates, (a different stimulus set in each block) with features differing on five cues. The same 12 pictures were presented 4 times for 5 seconds to each participant. After the learning phase, participants were presented again with all 12 old pictures as well as 20 new pictures in the testing phase. In each trial of the testing phase, participants had to indicate if the picture was an old picture or a new one.  These two phases were repeated for each of the five stimulus sets. 
 
## The Data set 
 
The data I will use contains the following variables:
 
- `ID`: participants ID
- `block`: block number 
- `trial`: trial number
- `stimulus`: type of stimulus (house, alien, plant, plate, or bug)
- `response`: response given in this trial (old or new)
- `correctResponse`: correct response in this trial
- `correct`: is `true` if response is equal to `correctResponse`, if not is it `false`
 

|ID               | block| trial|stimulus |response |correctResponse |correct |
|:----------------|-----:|-----:|:--------|:--------|:---------------|:-------|
|hcpibo8hk7bz2itt |     1|     1|houses   |old      |new             |false   |
|hcpibo8hk7bz2itt |     1|     2|houses   |new      |new             |true    |
|hcpibo8hk7bz2itt |     1|     3|houses   |old      |new             |false   |
 
Based on these variables, I calculated the number of hits, false alarms, false negatives, and misses, as well as their corresponding rates for each person, in each block. In addition, I calculated $d'$  as $d' = z(HR) - z(FR)$ (Snodgrass &  Corwin, 1988; Stanislaw & Todorov,1999), where HR again is the hit rate and FR is the false-alarm rate.
 

{% highlight r %}
hits <- dataSDT %>%  
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
 
Which then gives us the following data structure:
 

| IDn|stimulus | block| n_old| n_new|  h| fa|  hit_rate| fa_rate|  z_h_rate| z_fa_rate|     dprime|
|---:|:--------|-----:|-----:|-----:|--:|--:|---------:|-------:|---------:|---------:|----------:|
|   1|aliens   |     5|    12|    20|  7| 11| 0.5833333|    0.55| 0.2104284| 0.1256613|  0.0847670|
|   1|bugs     |     2|    12|    20|  7| 13| 0.5833333|    0.65| 0.2104284| 0.3853205| -0.1748921|
|   1|houses   |     4|    12|    20|  7| 14| 0.5833333|    0.70| 0.2104284| 0.5244005| -0.3139721|
|   1|plants   |     3|    12|    20|  7| 12| 0.5833333|    0.60| 0.2104284| 0.2533471| -0.0429187|
|   1|plates   |     1|    12|    20| 11| 11| 0.9166667|    0.55| 1.3829941| 0.1256613|  1.2573328|
 
 
### SDT - Non Hierarchical - One Stimulus Set 
 
I will first build the non-hierarchical version of the SDT model as described in Chapter 11 on pages 158-159 in the fantastic book of Lee and Wagenmakers (2014), analyzing the data from only one stimulus set. 
 
 

{% highlight r %}
temp <- hits %>% filter(stimulus == "plants")
{% endhighlight %}
 
#### The Model 
 
Next, I will define the SDT model in STAN.  
 
##### The Data 
 
In the `data` block, I define the data we are using in the model, which is:
 
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
 
stan_data
{% endhighlight %}



{% highlight text %}
## $s
##  [1] 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12
## [30] 12 12 12 12 12 12 12 12 12 12 12
## 
## $n
##  [1] 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20
## [30] 20 20 20 20 20 20 20 20 20 20 20
## 
## $h
##  [1]  7  9  7 10 11 12 10  7  6  7  7  8  8  7  6 12  6 11  9 12  7  6 10  8  8  6 11  5  7
## [30]  8  9 10  6  9  6  6  7  5  6  9
## 
## $fa
##  [1] 12 10 16 10 11 20 12 10 13  7 11 12 11 10  6 20 11 11 14 11 10 14 14 17 10 15  9 11 11
## [30] 13 16 10 12 13 16 11  8 11 11 13
## 
## $p
## [1] 40
{% endhighlight %}
 
##### The Parameters 
 
In the `parameters` block, I define the two main parameters we have in our model and which we are interested in, *d* and *c*. These two parameters are then transformed in to the *hit rate* and the *false-alarm rate* in the `transformed parameters` block. In the `model` block, I then define the binomial likelihood function connecting our parameters and our data. This can be written as: 
 
 
$$
\begin{aligned}
HR &= \phi(0.5 \times d - c) \\[.5em]
FR &= \phi(- 0.5 \times d - c) 
\\[1.5em]
h & \sim binomial(h,s)\\[.5em]
fa &\sim binomial(fa,n)
\end{aligned}
$$
 
 
 
 
In the `generated quantities` block I also include variables for later posterior predictive analysis, capturing the predictions of `h` and `fa` based on the current parameter values of each step of the MCMC-chains. 
 
##### The Priors 
 
The priors for both *d* and *c* are normal distributions with $\mu = 0$ and $\sigma = 1$
 
$$
\begin{aligned}
d &\sim Normal(0,1)\\[.5em]
c &\sim Normal(0,1)
\end{aligned}
$$
 
 
which corresponds to a uniform distribution after transforming them into hit and false-alarm rates.
 
 
<img src="/assets/img/2020-04-28-SDT1-nonhierarchical.Rmd/unnamed-chunk-4-1.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />
 
##### The STAN-Model 
 
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
 
#### Run the Model 
 
Next, we sample from the model using 10.000 iterations, with  a rather small warm-up of 2000 iterations and thinning = 4. 
 

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
 
# save results 
posterior_SDT_oS_nH <- rstan::extract(fit1, permuted = FALSE)
summary_oS_nH       <- summary(fit1) %>% as.data.frame()
{% endhighlight %}
 
#### Inspect MCMC, Rhat, ESS
 
First, we can look at some MCMC-Traces for some of the parameters and persons. 
 

{% highlight r %}
 mcmc_trace(posterior_SDT_oS_nH,
            pars = vars("d[1]":"d[3]",
                        "c[1]":"c[3]",
                        "hit_rate[1]":"hit_rate[3]",
                        "fa_rate[1]" :"fa_rate[3]"),
            facet_args = list(nrow = 4, labeller = label_parsed))
{% endhighlight %}

<img src="/assets/img/2020-04-28-SDT1-nonhierarchical.Rmd/unnamed-chunk-7-1.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" style="display: block; margin: auto;" />
 
So far so good, this looks exactly as you would like it, some nice hairy caterpillars. We can also plot the distributions of $\hat{R}$  and the *effective sample size*:
 

{% highlight r %}
p1 <- ggplot(summary_oS_nH,aes(x =summary.n_eff))+
        geom_histogram(bins = 40, color = "black", fill = "skyblue2")+
        theme_bw() +
        labs(x = "Effective Sample Size",
             y = "Count")
 
p2 <- ggplot(summary_oS_nH,aes(x =summary.Rhat))+
        geom_histogram(bins = 40, color = "black", fill = "tomato2")+
        theme_bw() +
        labs(x = "Rhat",
             y = "Count")
 
p1+p2 #patchwork package
{% endhighlight %}

<img src="/assets/img/2020-04-28-SDT1-nonhierarchical.Rmd/unnamed-chunk-8-1.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" style="display: block; margin: auto;" />
 
This looks also fine. The effective sample sizes are always > 400, which is often recommended, and are near the total sample size of  4000. Also $\hat{R}$ is always very close to 1.
 
#### Posterior Summaries & Densities
 
Next, we can look at the summary statistics and plots of the posterior density distributions. So lets make a convenient tidy data.frame for plotting and for the posterior summary statistics. 
 

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
 
mcmc_SDT_oS_nH <-  rstan::extract(fit1,pars="lp__",include=FALSE) %>%
                      bind_rows() %>% 
                      mutate(ID     = rep(1:40,each = 4000),
                             dprime = rep(temp$dprime, each = 4000)) %>% 
                      group_by(ID) %>% 
                      mutate(mPP_d  = median(d),
                             mPP_c  = median(c))
{% endhighlight %}
 
 
##### d 
 
Since *d* is the parameter we are most interested in, lets start with this one. Below you see a forest plot showing the median, 50%, and 96% credible intervals of the posterior distributions, as well as uni-variate marginal posterior distributions for some participants, showing the median posterior estimate (red dashed line) as well as the analytically calculated *d* (black dashed line).
 

{% highlight r %}
mcmc_intervals(posterior_SDT_oS_nH,
               pars = vars(starts_with("d")),
               point_est = "median",
               prob = 0.5,
               prob_outer = 0.96)
{% endhighlight %}

<img src="/assets/img/2020-04-28-SDT1-nonhierarchical.Rmd/unnamed-chunk-10-1.png" title="plot of chunk unnamed-chunk-10" alt="plot of chunk unnamed-chunk-10" style="display: block; margin: auto;" />
 
 

{% highlight r %}
mcmc_SDT_oS_nH %>%
  filter(ID %in% 1:6) %>% 
  ggplot(., aes(x = d))+
    geom_density(color = "black",fill = "skyblue2",alpha=0.5) +
    geom_vline(aes(xintercept = mPP_d), color = "red", lty = "dashed",lwd = 1) +
    geom_vline(aes(xintercept = dprime), color = "black", lty = "dashed") +
    facet_wrap(.~ID,scales="free") + 
    theme_bw()
{% endhighlight %}

<img src="/assets/img/2020-04-28-SDT1-nonhierarchical.Rmd/unnamed-chunk-11-1.png" title="plot of chunk unnamed-chunk-11" alt="plot of chunk unnamed-chunk-11" style="display: block; margin: auto;" />
 
We can see that the *d* values are rather small and very close to 0 for most people. This indicates that participants had a hard time differentiating old from new stimuli. I had hoped for larger values, as it is necessary that people are able to discriminate between stimuli and their features rather well for the multiple-cue judgment experiment I wanted to conduct with these stimuli. Also, from the plots it is evident that the median of the posterior distribution is very close to the analytically calculated *d'* value.  
 
 
##### c
 
We can also look at the *c* values. Since we have more noise trials (new stimuli) than signal trials (old stimuli) in our testing phase (20 vs. 12), and participant were told this information, I expected to find slightly positive values of *c*. However, as apparent from the summary statistics and the plots, participants had more negative values of *c*, indicating a bias for the "old"-response.
 

{% highlight r %}
mcmc_intervals(posterior_SDT_oS_nH,
               pars = vars(starts_with("c")),
               point_est = "median",
               prob = 0.5,
               prob_outer = 0.96)
{% endhighlight %}

<img src="/assets/img/2020-04-28-SDT1-nonhierarchical.Rmd/unnamed-chunk-12-1.png" title="plot of chunk unnamed-chunk-12" alt="plot of chunk unnamed-chunk-12" style="display: block; margin: auto;" />
 

{% highlight r %}
mcmc_SDT_oS_nH %>%
  filter(ID %in% 1:6) %>% 
  ggplot(., aes(x = c))+
    geom_density(color = "black",fill = "skyblue2",alpha=0.5) +
    geom_vline(aes(xintercept = mPP_c), color = "red", lty = "dashed",lwd = 1) +
    facet_wrap(.~ID,scales="free") + 
    theme_bw()
{% endhighlight %}

<img src="/assets/img/2020-04-28-SDT1-nonhierarchical.Rmd/unnamed-chunk-13-1.png" title="plot of chunk unnamed-chunk-13" alt="plot of chunk unnamed-chunk-13" style="display: block; margin: auto;" />
 
### The posterior predictive values
 
We can also look at the posterior predictive values of `h` and `fa`.  For this, I first combine and then plot the actual observed values from our original data.frame with the MCMC estimates. 
 

{% highlight r %}
temp     <- hits  %>%  filter(stimulus == "plants") %>% ungroup() %>% select(.,IDn,h,fa)
postpred <- mcmc_SDT_oS_nH %>% 
              select(., ID, h_pred, fa_pred) %>% 
              left_join(.,temp,by = c("ID" = "IDn")) 
 
postpred[1:5,]
{% endhighlight %}



{% highlight text %}
## # A tibble: 5 x 5
## # Groups:   ID [1]
##      ID h_pred fa_pred     h    fa
##   <dbl>  <dbl>   <dbl> <dbl> <dbl>
## 1     1      8      11     7    12
## 2     1      7      13     7    12
## 3     1      9      14     7    12
## 4     1      4      14     7    12
## 5     1      9      11     7    12
{% endhighlight %}
 
 

{% highlight r %}
pp1 <- ggplot(filter(postpred,ID %in% 1:4),aes(x = h_pred)) +
        geom_histogram(bins=10,color = "black",fill = "skyblue2",alpha=0.5) +
        geom_vline(aes(xintercept =  h),color = "black",lwd = 1.5, lty = "dashed")+
        facet_wrap(.~ID,scales="free") + 
        theme_bw() 
 
pp2 <- ggplot(filter(postpred,ID %in% 1:4),aes(x = fa_pred)) +
        geom_histogram(bins=10,color = "black",fill = "tomato2",alpha=0.5) +
        geom_vline(aes(xintercept =  fa),color = "black",lwd = 1.5, lty = "dashed")+
        facet_wrap(.~ID,scales="free") + 
        theme_bw() 
 
pp1 + pp2 #patchwork package
{% endhighlight %}

<img src="/assets/img/2020-04-28-SDT1-nonhierarchical.Rmd/unnamed-chunk-15-1.png" title="plot of chunk unnamed-chunk-15" alt="plot of chunk unnamed-chunk-15" style="display: block; margin: auto;" />
 
This looks also fine so far. The posterior distribution is symmetrically distributed around the empirical values. In addition, let me also calculate how often the 95% credible interval of the predicted values contains the true values:
 

{% highlight r %}
temp1 <- hits  %>%  filter(stimulus == "plants") %>% ungroup() %>% select(.,ID=IDn,h,fa)
 
temp2 <- desc_oS_nH_tidy %>% 
          select(param,ID,HDI025,HDI975)  %>% 
          filter(param == "h_pred" | param == "fa_pred") %>% 
          pivot_wider(.,id_cols    = ID,
                        names_from = param,
                                          values_from = c(HDI025,HDI975))
 
left_join(temp1,temp2,by = "ID")  %>% 
  summarize(pp_h  = mean(h > HDI025_h_pred & h < HDI975_h_pred),
            pp_fa = mean(fa > HDI025_fa_pred & fa < HDI975_fa_pred)) %>% 
  as.data.frame() # better output in .md
{% endhighlight %}



{% highlight text %}
##    pp_h pp_fa
## 1 0.925  0.95
{% endhighlight %}
 
As evident from the plots and the 95% credible interval checks of the posterior predictives, our model is able to recover our data well.
 
---
 
### References 
 
- Lee, M. D., & Wagenmakers, E. J. (2014). *Bayesian cognitive modeling: A practical course.* Cambridge university press.
 
- Snodgrass, J. G., & Corwin, J. (1988). Perceptual Identification Thresholds for 150 Fragmented Pictures from the Snodgrass and Vanderwart Picture Set. *Perceptual and Motor Skills*, 67(1), 3â€“36. https://doi.org/10.2466/pms.1988.67.1.3
 
- Stanislaw, H., & Todorov, N. (1999). Calculation of signal detection theory measures. *Behavior Research Methods, Instruments, & Computers*, 31(1), 137â€“149. https://doi.org/10.3758/BF03207704
 
 
