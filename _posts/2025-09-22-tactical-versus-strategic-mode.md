---
layout: post  
comments: true  
tags:   
title: Strategic Versus Tactical Mode for AI/ML-enabled Solution Development Teams
---

Data science and AI engineering teams need to be able to operate in two distinct but complementary modes in order to design, build, and deploy AI/ML-enabled solutions effectively: strategic mode and tactical mode. Lately, I’ve been reflecting on how these modes show up, both in my own work as an individual contributor and in leading teams. My thinking has been inspired by Daniel Kahneman’s concepts of System 1 and System 2 thinking from his book *Thinking, Fast and Slow*.

In practice, I’ve noticed that teams often lean more heavily into one mode depending on the nature of the work and the stage of the product lifecycle. While it’s an oversimplification to say that one mode maps cleanly to one phase, the framework is useful for understanding what teams need and when they need it.

## Strategic Mode  
What does it mean for a data science or AI engineering team to be in strategic mode? It’s the expansive phase, when we step back and consider not just what we’re building, but how it fits into the broader ecosystem of team goals and business objectives.

In this mode, the focus is on questions like:

* How can we design for generalizability and modularity?
* Is there a platform perspective we can take?
* Could this solution serve multiple use cases beyond the immediate one?
* How do we build something today that makes future development, testing, and deployment easier?

Strategic mode is about positioning work in a way that maximizes long-term value—for both end users and the organization. It’s closer to Kahneman’s System 1: broad, intuitive, fast-moving insights about an expansive solution space.

Most teams naturally lean into strategic mode during the design and planning phase. This is where opportunities for “productivity multipliers” emerge. These are things like internal tooling, reusable modules, or even full-fledged platforms. Many of the most powerful enablers of future efficiency come from investing in strategic thinking up front.

To succeed in this mode, teams need:

* **Time** to explore the solution space and build prototypes.
* **Space** for brainstorming and collaboration.
* **Information** about user needs, technical constraints, and business priorities.

At this stage, open-ended questions are especially useful. Asking about modularity or generalizability can uncover opportunities that otherwise might be overlooked.

## Tactical Mode
Tactical mode is where the focus narrows. This is the detailed, methodical, “in the weeds” phase of work. Here, teams are less concerned with broad positioning and more with making decisions that enable progress day-to-day.

Typical questions in tactical mode include:

* Which vector database should we use for our RAG system?
* How do we design a clear audit trail for regulatory compliance?
* What’s the best framework for this specific use case?

This mode aligns more with Kahneman’s System 2: deliberate, systematic, and detail-heavy.

While strategic thinking dominates design and planning, tactical thinking takes center stage during building, testing, and deploying MVPs and beyond. It’s about speed, precision, and quality—navigating issues that inevitably arise when a prototype becomes a production system.

To thrive in tactical mode, teams benefit from structure. Ideally, design and planning outputs include:

* A clear high-level **architecture** (inputs, outputs, integrations, and regulatory considerations).
* A rough **sprint plan** and buffer time to account for uncertainty.
* A **target date** for the first release.

This structure channels effort toward execution—and execution is inherently tactical.

## Reality
Of course, saying that planning is strategic and building is tactical is too neat. In reality, both modes are in play throughout the lifecycle. During planning, teams still make tactical choices about fit-for-purpose technologies. During execution, teams often surface strategic insights—like discovering a reusable solution or identifying opportunities for internal tooling.

The balance shifts over time: strategy dominates early, tactics dominate during execution. But keeping both modes in mind helps leaders create the right conditions for their teams—whether that means carving out space for big-picture thinking or providing structure to support focused execution.