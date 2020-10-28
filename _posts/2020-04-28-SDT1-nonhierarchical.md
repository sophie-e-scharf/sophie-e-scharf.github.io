---
layout: single
title: "Signal Detection Model with STAN - Part 1"
date: 2020-04-28
tags:
  - cognitive modeling
  - STAN
  - R
  
# status: process
published: true
status: publish
# published: false
permalink: /posts/2020/04/SDT-1/
---
 

 
 
 
## About this blog post
 
This blog post is supposed to be the first in a series of blog post comparing the *signal detection theory * model (SDT) with a *two-high-threshold* model (2HTM) of recognition. I wrote this blog post for three reasons. First, I wanted to learn more about the modeling of recognition and memory data. Second, I wanted to learn more about STAN, since I mostly use JAGS in my own research. Finally, I also wanted to practice writing, since I take forever when writing my own articles. So the reasons for this blog are rather selfish. However, if anyone ever finds this blog post and finds it helpful, that would be even better!
 
 
In this blog post you are going to read, I will use a *non-hierarchical SDT* model to investigate the data from a recogntion experiment. In following blog posts, I will extend this model to account for differences between individuals as well as differences between stimulus sets. Then I will model the same data with a 2HTM and finally compare both models with each other. 
 
## Setup 
 
At the beginning, I load the packages I need for this analysis. However, before we start with the actual analysis and modeling, I first want to give a short introduction into signal detection theory.
 










































