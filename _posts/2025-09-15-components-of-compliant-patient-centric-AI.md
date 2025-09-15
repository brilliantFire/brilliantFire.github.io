---
layout: post  
comments: true  
tags:   
title: Components of Compliant, Patient-Centric AI Strategy
---

Patient safety and satisfaction should be at the heart of AI intiatives in the clinical space. At a bare minimum, AI strategies need to include consideration for the regulatory requirements invloved in initiatives that directly or indirectly affect patient safety, care, or satisfaction. To be truly patient-centric, more is required to ensure that we are measuring outcomes to correctly assess the impact of an AI solution.

In this post, we explore possible components of a compliant, patient-centric AI solution strategy. We'll break them down into things that should be done or considered at the beginning of AI intiatives, during the development of the solution, and in post-production.

## Before development begins

### Planning for the *entire* project
Planning is the one of the biggest, most important components of a compliant, patient-centric AI solution strategy. In fact, we could easily break planning down further by what the planning is for (development, testing, or end-to-end process validation to meet regulations). All AI initiatives do some form of planning, obviously, but many times these efforts stop short of the entire product lifecycle. This is especially true for what comes after solutions are put into production. Things like incident management and monitoring processes should be worked out ahead of time so that solutions can hit the ground running. At the planning stage, any validation plans should be made to the extent that they can. We should also have identified all the appropriate stakeholders and teams required to carry out the initiative (engineering, product, and business). There's a lot to consider beyond the technical design of the solution.  

### Metrics that measure what matters most
Metrics help us to understand whether or not the AI initiative is successful. A big part of that is patient impact. Metrics and their measurement should be a topic of planning sessions and measuring needs to start well beofre development does. Presumably, there exists some intake process by which AI initiatives are vetted and prioritized. Out of this process, there should be at least some intial measurements to guide planning. If not, these metrics need to start being collected during planning. What sorts of metrics should be collected? These should heavily skew towards the patient-centric. Depending on the use case, things like time-to-appointment/treatment, wait time in clinic, or patient-reported outcome measures should be measured to understand the process the AI initiative will be a part of. It is not enough to meaure outcomes, though. Patient satisfaction is a key metric that is often overlooked. Efforts should be made to measure this as well.

### Internal and external review
Once plans have crystalized for the development, testing, and validation of an AI solution, these plans should be reviewed by people internal and external to the team. Internally, there should be alignment among all the groups involved in the initiative, including various engineering groups, product, and business stakeholders. Every patient-impacting initiative should consider an ethics review by a committee consisting of end users, patients, and possibly other groups like legal. To ensure compliance with regulations, solutions should also conduct thorough quality assurance and regulatory reviews.

## During development

### Consistent, frequent communication
Communication does not end once plans are finalized and reviewed. The technical team doing the active development will ideally be tightly iterating on solutions with product teams and the business. At Moderna, our most successful projects had regular meetings among a core group of technical, product, and business stakeholders. The frequency of the meetings ranged from weekly during planning stages to 2-3 times per week during execution to daily during post-production hypercare. Keeping in close communication with everybody meant that there were no surprises during development. We were all on the same page as to where we were in the process and how well things were moving.  

### Comprehensive testing
A comprehensive testing strategy should also be developed in the planning stage. This strategy starts with unit, integration, and end-to-end testing of the technical solution. Also during the early development stage, any machine learning or AI components should be backtested by stratifying predictions according to sensitive variables, such as race, gender, or age. Before going to production, the solution will need to have user acceptance testing completed. Finally, any official system validation process, like for GxP processes in pharma, needs to be completed. 

## Post-production

### Monitoring for bias and data drift
If our planning was as thorough as it needed to be, it should be clear what needs to be monitored and why. For patient-impacting AI initiatives, bias monitoring is a must. This can be achieved through periodic backtesting against sensistive groups, or in an automated fashion by [calculating a metric such as disparate impact](https://medium.com/ibm-data-ai/fairness-explained-definitions-and-metrics-9690f8e0a4ea). Bias monitoring is often spoken of together with drift monitoring which monitors for changes in the distribution of data between what models are trained on and the production data they score. All metrics should be easily reportable. I deally, there is a dashboard that visualizes  updates to these metrics as they are measured. This is one area where platform thinking can really make a difference. Teams might consider developing generalized tooling to help them calculate metrics and monitor models. This tooling could be used throughout the development and post-production periods and used across different initiatives.

### How to handle incidents
Again, plans for incident management once the system goes into production should be made early and adjusted as needed throughout development. Incident management needs to be taken very seriously in regulated systems. Incidents and their resolutions need to be closely tracked. Ideally, incident management is handled through a system that end users and other key stakeholders can easily access so that they can report incidents in a timely fashion and monitor progress in solving them.

### User and/or patient satisfaction
To be patient-centric, we must monitor patient satisfaction for processes that include AI or ML components. This might mean looking beyond the end user if the end user is internal to the organization providing the patient service. For example, if we develop a solution to assist hospital staff with scheduling medical procedures, we should ask both hospital staff and patients that have been scheduled to score the process. Creatively finding ways to measure patient satisfaction for AI initiatives that are indirectly affecting it could reveal shortfalls that might otherwise be overlooked.

## Bringing it all together
A patient-centric AI strategy requires attention to the full lifecycle of an initiative: before development, during development, and after a solution is in production. Success depends on careful planning that prioritizes compliance and patient impact, defining metrics that measure safety and satisfaction, and reviewing plans with the right internal and external stakeholders. During development, consistent communication and thorough testing help ensure both quality and fairness. Once in production, monitoring, incident management, and collecting feedback on patient satisfaction keep solutions reliable and impactful. Taken together, these practices build AI initiatives that not only meet regulatory requirements but also improve the patient experience in meaningful ways.
