---
layout: post
comments: true
tags: agile change-management data-science-management data-sciences-operations leadership
title: Update, Changes, and Adopting Agile  
---

Welcome back! It's been 2,151 days since my last post. A lot has happened. Lucky for you, I won't talk about the vast majority of it here (who wants to relive the pandemic, amirite?) but I figure I'll catch you up a little on where I'm at professionally, what it's like working at a company that's undergoing a lot of growth and changes ([Moderna](https://www.modernatx.com/en-US)), and a direct consequence of that growth for my team - moving to doing agile data science.  Let's get started!  

### What's up?  
I'm currently a Senior Principal Data Scientist at Moderna, a biotech company specializing in mRNA medicines for infectious disease, oncology, and rare disease. I lead a team of four data scientists. Together, we comprise Data Science & Artificial Intelligence for Development (DSAI4D). We work on AI-enabled software for Moderna's 1,200-person-strong Development organization. If you're familiar with the [drug development lifecycle](https://www.fda.gov/patients/learn-about-drug-and-device-approvals/drug-development-process), Development covers everything after pre-clinical studies but before commercial marketing. Among other things, we help the teams at Moderna who:

* Write our Investigational New Drug and Biologics License Applications  
* Design and conduct our clinical trials  
* Manage and analyze the data collected during clinical trials  
* Follow up on adverse events both during clinical trials and after one of our drugs goes to market  

Our team is part of a larger Digital for Development organization consisting of technical program managers, product managers, data engineers, and software engineers. We roll up to the Chief Information Officer, currently Brad Miller.  

### Ch-ch-ch-Changes  
I have been at Moderna for a little less than three years and during that time the company has grown significantly. A consequence of this growth has been constant change as our organization morphs from a small biotech company with less than 1,000 people in 2019 to a full-fledged medicine company currently employing over 5,000. Part of this dramatic growth has been a re-committment to one of Moderna's original ambitions - to be the ["digital biotech company"](https://www.modernatx.com/en-US/media-center/all-media/blogs/building-the-first-digital-biotech). Functionally, I like to think of us as striving to have a tiny tech company inside a biotech company. Moderna was founded with the idea of "digitizing everywhere". This includes all the necessary elements to build first-class AI-infused software for the sole purpose of helping our Development organization get powerful mRNA medicines to people faster.  

### How I Learned to Relax and Adopt Agile  
Let me just begin by saying that there's no one right way to "do" data science on an organiational level. Each organization will have to assess their own specific needs. In some organizations, data science takes on more of a research quality, where the goal is to support the business with sophisticated analyses that are consumed in dashboards, chatbots, or small apps. For small organizations where data scientists might where many hats, this might be the best way to maximize the value data scientists provide. Initially starting the DSAI4D team at Moderna, I took this approach to distributing workload across the team. I ran the team much like one would run a reasearch lab - single scientists going deep on one or two projects without much collaboration with other data scientists.

For larger organizations, however, having "full-stack" data scientists own 100% of the project lifecycle is an inefficient way to utilize a data scientist's skills given that, with adequate planning, some of the work can be distributed to data and software engineers, allowing scientists to focus on model building and evaluation. Furthermore, in situations like ours at Moderna where the digital organization needs to support a large number of stakeholders with a small number of developers, it is simply untenable to have one scientist tied up on a single product. Finally, single-threading products on the team is unstable and lacks the redundancy required to keep product development moving when their primary developers are away.  

To take advantage of all these things, and to really create that "tech company inside a biotech company", our data science team needed to adopt [agile methods of working](https://www.atlassian.com/agile). Currently, we are working in a "scrumban" framework. This combination of [scrum](https://www.atlassian.com/agile/scrum) and [kanban](https://www.atlassian.com/agile/kanban) works well for teams like ours that are:

1. Just making the transition to agile, and;  
2. Having a lot of uncertainty early on in product development.   

Here are the key takeaways I've learned so far in the process of moving from "research lab" to "agile development team":  

#### Relax About Uncertainty and Time  
One of my biggest early complaints about scrum and data science - and one I think shared by a lot of practitioners - was how to deal with the greater level of uncertainty that we often have with data science projects or features. The trick here for me was to RELAX. It turns out agile and even scrum are flexible enough to handle even data science-levels of uncertainty. Reframing answers to questions or the results of experiments as deliverables can help track progress towards design goals. "Spike" issues can be used to track this sort of work. Don't worry about 2 week sprint lengths, or trying to time-box work early on. The work takes as long as it takes. Agile, scrum, kanban. They are all just ways of tracking that work.  

#### Get Used to Tracking Work  
Agile - with its issues and t-shirt-sizing and other informational elements - can sometimes appear to add an additional documentation burden for scientists. I've found that this gets easier and lighter with experience. Work starts with an issue - this could be a ticket in Jira or a row in a Google Sheet - and ends in the issue with documentation of results or work completion. Writing "good" issues with detailed descriptions and acceptance criteria can make it easier later when it comes time to write up results. Over time, I've been able to add working within issues to my regular data science workflow. Now it's just as natural as opening VS Code.

#### Really Lean-In to the Technology  
There are a lot of platforms for tracking agile work. One of the most popular is Jira. We've had a lot of success leveraging this technology to track work and do reporting on a large organizational level. Whatever it is you are using (even Google Sheets) could have templates or other features to help track work and report on it. Finding out and applying features in tracking technology can lower the burden of documentation and make transitioning to an agile way of working easier.

Remembering that agile is time-flexible, getting used to tracking work in issues, and adopting technology to make tracking work easier have all helped DSAI4D so far on its agile journey. 

Stay tuned for more thoughts. I promise it won't be another 6 years until the next one ;)  

