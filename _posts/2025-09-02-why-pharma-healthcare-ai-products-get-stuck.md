---
layout: post  
comments: true  
tags:   
title: Why AI Products Get Stuck in Pilot Stage in Pharma and Healthcare
---

A recent study by a group a MIT revealed that a staggering [95% of AI pilots fail](https://fortune.com/2025/08/21/an-mit-report-that-95-of-ai-pilots-fail-spooked-investors-but-the-reason-why-those-pilots-failed-is-what-should-make-the-c-suite-anxious/) to deliver any savings or profits. Among the chief reasons underlying this observation are a learning gap where people and organizations fail to understand how best to use AI technologies and a lack of talent to build in-house systems or properly configure and run off-the-shelf ones.  

While these problems could be endemic to any organization in any industry, pharma and healthcare present special cases where the success of an AI product might have even less of a chance of really changing things for workers and unlocking value. Here we'll discuss a few of the special circumstances surrounding AI initiatives in the pharma and healthcare settings.

## Why Pilots are Easy (and Attractive)
Let's face it. Anybody can build a prototype. Especially now in the era of vibe-coding and low-code tools for making and deploying things like agentic workflows ([n8n](https://n8n.io/ai/) is a good example), coming up with a prototype relatively quickly has become at least 10x easier. Pilots are low-risk and low-cost, making them attractive to executives looking to "do AI". All too often, though, these pilots are not fully thought-through and quickly run into issues when trying to scale. Leadership is often satisfied with plans for the pilot without having clear plans for scaling.

## Challenges in Scaling in Pharma and Healthcare
The first challenges encountered when trying to scale up AI pilots are often a lack of good user requirements and the absence of a solid business case outlining ROI. Because of how easy it is to build one, a pilot is often constructed before the time is taken to really understand the problem. Upstream and downstream systems and business processes that the AI solution needs to integrate with are often overlooked. Another scale-up challenge is infrastructure. In the past, I have encountered scale-up issues because our in-house infrastructure was being constructed alongside our use cases. While this had the advantage of ensuring that our infrastructure met our needs precisely, it did sometimes slow down our ability to scale up.

While these techinical reasons can impede the progress toward production, so can cultural and regulatory reasons. Users such as clinicians, pharma safety monitors, and regulators are often (and with good reason) wary of black-box models. Special attention needs to be paid to the concerns of these individuals since the success of an AI initiatives critically depends on their expertise and insight. Often these conincide with another frequent roadblock: regulations. Pilots are not available for use in GxP processes or are considered standard-of-care because they don't go through the thorough processes that would make them useable in those contexts.

## How to Move from Pilot to Impact
So how do we move from a successful pilot to something that really starts generating value? A lot of that work should happen before the pilot is even started. First, crucial alignment between clinicians or end users on every detail of the proposed solution needs to be acquired. Here's where tools like [Lovable](https://lovable.dev/) can help give end users a glimpse of what the actual solution might look and feel like. User requirements should be detailed and documented in a system like Jira and they should include such critical things as any potential validation that would be required as a result of GxP or standard-of-care regulations. These early conversations begin to build cross-functional ownership: end users should feel like they are making contributions to the development of the system. In a similar sense, governance and monitoring should be discussed ahead of time and built in to the product development plan.

In addition to alignment and cross-functional ownership, a phased roadmap with clear milestones is necessary to tactically lay out how the product should be developed, tested, and rolled out to users. This should include early iteration on complete end-to-end processes together with end users to refine the product prior to a broader launch.

## A Playbook for Leaaders in Pharma and Healthcare
Taken together, these points boil down to a simple overall playbook to ensure AI products don't get stuck in the pilot phase.

1. Start with the problem (not the technology).  
2. Define ROI from the start (rough but realistic estimates are fine).  
3. Build infrastructure and compliance foundations (to the extent that we can predict what we need).  
4. Engage stakeholders across functions (and do so early and often).  
5. Scale only after *measurable* pilot success (design the pilot so it solves the problem as end-to-end as possible)

In the end, scaling AI in healthcare and pharma is about bridging the vision with execution. As we are seeing more and more in industries across the board, successful organizations are those that see AI as a key enabled of their business processes, not just a side project.