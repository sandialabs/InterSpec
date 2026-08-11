# InterSpec LLM Knowledge Directory

This directory holds local "skills" for the DeepResearch sub-agent.

For InterSpec today, these skills are primarily curated knowledge packets (not full agent capabilities). The goal is to inject relevant, trusted context into the LLM when DeepResearch is asked domain questions.

## Why this exists

- Lets you add organization- or workflow-specific knowledge without recompiling InterSpec.
- Helps DeepResearch use curated information instead of relying only on model memory.
- Keeps local expert notes and procedures

## Where to put skills

InterSpec looks for `llm_knowledge` in two places:

- App/static data directory (this repository's `data/llm_knowledge`).
- User writable data directory (recommended for your own skills).

To find your writable data directory in InterSpec:

1. Open **About InterSpec**.
2. Go to the **Data** tab.
3. Look in the **user data** section.

Typical writable data locations:

- Windows: `C:\Users\<username>\AppData\Roaming\InterSpec`
- macOS: `/Users/<username>/Library/Containers/gov.sandia.macOS.InterSpec/Data/Library/Application Support/sandia.InterSpec`

Create a `llm_knowledge` directory under that directory if it does not already exist.

## How to add a skill (knowledge packet)

Create a subdirectory per topic under `llm_knowledge`, and put a `SKILL.md` file inside it.

Example:

```text
llm_knowledge/
  477keV_feature/
    SKILL.md
  neutron_transport_notes/
    SKILL.md
```

`SKILL.md` can use Agent Skills style frontmatter:

```yaml
---
name: neutron-transport-notes
description: Notes on moderation, capture, and interpretation pitfalls.
---
```

Everything after frontmatter is the knowledge body that may be loaded into DeepResearch context.

For the broader format and ecosystem, see:

- [https://agentskills.io/home](https://agentskills.io/home)

## Current implementation limits

Current behavior is intentionally simple:

- Only `SKILL.md` is used.
- Only `name` and `description` frontmatter fields are read (best effort).
- DeepResearch can load the `SKILL.md` body as context; other Agent Skills features are not implemented yet.
- `scripts/`, `references/`, and `assets/` are not executed/loaded as full skills support.
- Skills are exposed only through the DeepResearch path (not all agents directly).

In short: treat these as curated knowledge files for context injection, not as fully executable agent skills.

