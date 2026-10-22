# InterSpec LLM Knowledge Directory

This directory holds local knowledge "skills" for the KnowledgeBase sub-agent of the LLM Assistant.

For InterSpec today, these skills are primarily curated knowledge packets (not full agent capabilities). The goal is to inject relevant, trusted context into the LLM when the KnowledgeBase agent is asked domain questions.

## Why this exists

- Lets you add organization- or workflow-specific knowledge without recompiling InterSpec.
- Helps the LLM use curated information instead of relying only on model memory.
- Keeps local expert notes and procedures

## Where to put skills

InterSpec looks for skills in these directories, in increasing order of precedence:

1. `<static data dir>/llm_knowledge` (this repository's `data/llm_knowledge`)
2. `<static data dir>/gamma_spec_skills/llm_knowledge`
3. `<writable data dir>/llm_knowledge` (recommended for your own skills)
4. `<writable data dir>/gamma_spec_skills/llm_knowledge` (e.g., a `gamma_spec_skills` knowledge repository
   cloned directly into the writable data directory)
5. `<LLM config dir>/llm_knowledge`, and `<LLM config dir>/gamma_spec_skills/llm_knowledge`, if InterSpec
   was started with `--llm-config-dir` (as the MCP agent-server container is)

A skill with the same `name` as one in an earlier directory replaces it (including its `references/`).

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
  neutron-notes/
    SKILL.md
  detector-models/
    SKILL.md
    references/
      identifinder-family.md
```

`SKILL.md` can use Agent Skills style frontmatter:

```yaml
---
name: neutron-transport-notes
description: Notes on moderation, capture, and interpretation pitfalls.
---
```

Everything after frontmatter is the knowledge body that may be loaded into the KnowledgeBase agent's context.

For the broader format and ecosystem, see:

- [https://agentskills.io/home](https://agentskills.io/home)

## References (progressive disclosure)

A skill may have a `references/` sub-directory containing additional `.md` files:

```text
llm_knowledge/
  detector-models/
    SKILL.md
    references/
      identifinder-family.md
      detective-family.md
```

The skill's tool returns the `SKILL.md` body plus a list of `available_references` (name and one-line summary); the LLM can then call the tool again with `{"reference": "identifinder-family"}` to load only the detail it needs.  The summary is taken from an optional frontmatter `description:` field in the reference file, or else its first Markdown heading.

Keep `SKILL.md` short (a few KB): when to use the skill, the most important decision rules, and an index of what each reference covers.  Keep each reference self-contained and reasonably small (under ~10 KB); this matters for smaller LLMs.

## Benchmark question sets (optional)

A skill may carry its own LLM benchmark, so it can be evaluated wherever it is installed:

```text
llm_knowledge/
  detector-models/
    SKILL.md
    references/
    test_datasets/
      detector-models/
        detector_models_llm_benchmark.xml
        spectra/...            (spectrum files the XML references, by relative path)
```

The LLM Assistant's "Benchmarks / Quizzes" menu lists every `*_llm_benchmark.xml` found under the knowledge
directories (as well as under `test_datasets/` next to the data directories).  `SpectrumFile` paths in the XML
are relative to the XML file.  Do not put a `SKILL.md` under `test_datasets/` (it would be loaded as a skill).

## Current implementation limits

- Only `name` and `description` frontmatter fields are read (best effort).
- Only `SKILL.md` and `references/*.md` (not nested) are used; `scripts/` and `assets/` are not executed/loaded.
- Skills are exposed only through the KnowledgeBase sub-agent (`invoke_KnowledgeBase`), not to other agents directly.  That agent (and its `invoke_KnowledgeBase` tool) is only offered when at least one skill is found, or a `KnowledgeBaseRagUrl` is configured; the calling agents see only the skill `name`s (so make the name descriptive, e.g. `detector-models`), while the KnowledgeBase agent sees each skill's `description` when choosing which skill to load - so write a keyword-rich `description` (≤ ~400 characters).

In short: treat these as curated knowledge files for context injection, not as fully executable agent skills.
