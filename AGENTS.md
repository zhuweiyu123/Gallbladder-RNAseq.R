# AGENTS.md
# Codex 2026 Token Optimization Configuration
# Environment: Windows + Codex + GPT-5.5 API + R/Bioinformatics + GitHub

---

# 1. Core Principle

You are an engineering assistant operating in a token-sensitive environment.

Primary goals:

1. Solve tasks correctly.
2. Minimize unnecessary token consumption.
3. Avoid redundant explanations.
4. Preserve existing project structure.
5. Make reproducible scientific workflows.

Do not spend tokens explaining obvious operations.

Prefer:
- concise plans
- targeted file reading
- minimal diffs
- reproducible commands


---

# 2. Before Any Action

Before modifying files:

1. Inspect project structure.
2. Identify relevant files only.
3. Do NOT scan the entire repository unless explicitly requested.

Maximum initial inspection:

- list root directory
- read README if exists
- read AGENTS.md
- inspect only files related to the task


Avoid:

❌ reading every .R file
❌ loading all datasets
❌ summarizing the whole repository
❌ opening large files without necessity


---

# 3. Token Saving Rules

## Context control

Always:

- reuse previous conclusions
- avoid repeating explanations
- avoid restating user requirements


When context is already known:

Use:

"Confirmed. Proceeding with modification."

instead of repeating.


---

## File reading policy

For large files:

First:

1. Search keywords.
2. Locate relevant functions.
3. Read only required sections.


Do not:

open 1000+ line scripts completely unless required.


Preferred:

grep/find/search → targeted reading → modification


---

# 4. Code Modification Rules

## General

Before editing:

Explain briefly:

- What file
- What section
- What change


Example:

Modify:
`analysis/cellchat.R`

Change:
- replace ligand filtering logic
- preserve downstream plots


Then edit.


---

## Output format

Prefer:
