---
name: math
description: Update tex/CLAUDE.tex with mathematical insights from the current session
allowed-tools: Read, Write, Edit, Glob, Grep
---

# /math — Record mathematical progress

You are updating `tex/CLAUDE.tex`, a living LaTeX document that tracks the state of the research.
Think of it as the mathematical counterpart to `CLAUDE.md`: where CLAUDE.md tells Claude how to
*use* the codebase, CLAUDE.tex tells a *mathematician* what we know, what we've tried, and what's open.

## What to do

1. **Read** the current `tex/CLAUDE.tex` (if it exists).
2. **Review** the conversation so far. Identify mathematical content worth recording:
   - New results, conjectures, or observations (with proof sketches if discussed)
   - Computations or experiments and their outcomes
   - Dead ends and why they failed (these are valuable — record the *reason*, not just the failure)
   - Open questions or next steps that emerged
3. **Address CAPS flags.** Scan the document for comments or notes written in CAPITAL LETTERS (e.g., `% TODO CHECK THIS`, `% VERIFY`, `% IS THIS RIGHT?`). These are flags from the author requesting attention. For each one:
   - If the current session's work resolves the issue, fix the text and remove the flag.
   - If you can answer or clarify the question from context, do so and remove the flag.
   - If the flag is still unresolved, leave it in place but note any partial progress.
   - Summarise which flags you addressed (and which remain) in your report to the user.
4. **Update** `tex/CLAUDE.tex`:
   - **Add** new material in the appropriate section.
   - **Delete** or move to "Dead ends" anything the session showed doesn't work.
   - **Rearrange** freely — group related ideas together, promote important results, demote speculative remarks.
   - **Refine** existing entries if the session clarified or corrected them.
5. **Show** the user what you changed (brief summary).

## Style guide

- Write for a mathematician (the author). Assume familiarity with the notation in the research.
- Use standard LaTeX math: `$...$`, `\[...\]`, `align*`, `\begin{theorem}`, etc.
- Keep entries concise — a few lines each. This is a *log*, not a paper. Think: annotated bullet points with math.
- Date-stamp new entries with `% YYYY-MM-DD` comments so the evolution is traceable.
- Use `\section` for major categories, `\subsection` for themes within them.
- Cross-reference scripts and tex files by name when relevant (e.g., "see \texttt{scripts/verify\_conjecture.py}").

## Document structure

If creating `tex/CLAUDE.tex` from scratch, use this skeleton:

```latex
\documentclass[11pt]{article}
\usepackage{amsmath,amssymb,amsthm}
\usepackage[margin=1in]{geometry}
\usepackage{hyperref}

\newtheorem{theorem}{Theorem}
\newtheorem{conjecture}[theorem]{Conjecture}
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{observation}[theorem]{Observation}
\theoremstyle{definition}
\newtheorem{definition}[theorem]{Definition}
\theoremstyle{remark}
\newtheorem{remark}[theorem]{Remark}
\newtheorem{question}[theorem]{Question}

% Notation shortcuts (match the paper)
\newcommand{\qbin}[2]{\genfrac{[}{]}{0pt}{}{#1}{#2}_{\!q}}
\newcommand{\qpoch}[2]{(#1;q)_{#2}}
\newcommand{\GK}{\mathrm{GK}}
\newcommand{\CPP}{\mathrm{CPP}}
\newcommand{\QSym}{\mathrm{QSym}}

\title{Research Log: Cylindric Partitions and $A_2$ Andrews--Gordon Identities}
\author{Robin Langer}
\date{Last updated: \today}

\begin{document}
\maketitle
\tableofcontents

\section{Key results}
% Established facts, verified computations

\section{Conjectures and open questions}
% Active conjectures, things to try next

\section{Observations and patterns}
% Interesting but unproven patterns, experimental data

\section{Dead ends}
% Approaches that were tried and failed, with reasons

\end{document}
```

Adapt the sections as the research evolves. Add or merge sections as needed.

## Important

- This file is *yours to edit freely*. Delete stale content, correct errors, rewrite for clarity.
- If the session had no mathematical content worth recording, say so and don't make empty changes.
- If $ARGUMENTS is provided, treat it as guidance for what to focus on (e.g., `/math cyclic descent` means focus on that topic).
