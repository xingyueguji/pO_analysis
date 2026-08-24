---
name: combine
description: Use when working with CMS Combine (the HiggsAnalysis-CombinedLimit statistical analysis tool used in CMS searches and measurements). Covers datacards, physics models, limits, fits, running modes (AsymptoticLimits, FitDiagnostics, MultiDimFit, HybridNew, Significance, GoodnessOfFit, ChannelCompatibilityCheck), Combine error messages and warnings, statistical methodology, interpreting Combine output, and running Combine commands to reproduce, diagnose, or confirm results. Routes questions to the combine retrieval sources (docs / paper / code / forum) and runs Combine commands — in the shell when Combine is on PATH, otherwise via the remote combine-run server.
license: MIT
---

# Combine assistant

This skill covers two capabilities for CMS Combine
(HiggsAnalysis-CombinedLimit):

- **Answering questions** using the `combine` MCP server — read-only
  retrieval over four sources (`search_docs`, `fetch_doc`).
- **Running Combine** — either directly in the shell (when Combine is on
  the user's `PATH`) or, failing that, via the remote `combine-run` MCP
  server's `run_combine` tool.

Most requests are answered by retrieval alone. Reach for execution when
a user reports a command that misbehaves, asks you to run or try
something, or when reproducing a result would confirm a diagnosis.

---

# Part 1 — Answering questions (retrieval)

Use the `combine` MCP server. It exposes four complementary sources
through `search_docs` and `fetch_doc`. Route each question to the right
source(s), iterate intelligently, and answer with citations.

## The corpus

| Source ID | Covers |
|---|---|
| `combine-docs` | Official docs (MkDocs). How-to, CLI flags, tutorials, reference. |
| `combine-paper` | arXiv:2404.06614v2. Methodology, definitions, the "why". |
| `combine-code` | Source tree (v10.6.0). The implementation. |
| `combine-forum` | cms-talk Statistics category (2022→). Errors, workarounds, real-world edge cases. |
| `combine-hypernews` | Archived HyperNews forums (pre-2022). Same kind of Q&A, but older — **lower priority** than `combine-forum`; may reference outdated versions. |

## How to answer

### Step 1 — Pick the right source first

- **"How do I X?"** / **"What does flag --Y do?"** → start with
  `combine-docs`.
- **"Why does X work this way?"** / **"What is the formal definition of Y?"**
  → start with `combine-paper`.
- **"I'm getting this error"** / **"I see this warning"** → start with
  `combine-forum`; if it comes up empty, try `combine-hypernews` (the
  older archived forum — but treat version-sensitive answers with care,
  it predates 2022).
- **"What does the implementation actually do?"** / **"Is feature Z really
  there?"** → start with `combine-code`.

Sources are complementary, not redundant. If the first source comes up
empty, fall through to the next most likely one.

### Step 2 — Search and read scores intelligently

Call `search_docs(query, source)` to get up to 10 ranked hits.
Interpret the scores:

- **Top score is high (>15) with a clear gap to the runner-up:** the
  top hit is almost certainly the right one. Fetch it.
- **Top hits cluster tight (within ~3 points):** all of them are
  relevant. Skim 2–3 snippets to pick the best one before fetching.
- **All scores low (<8):** the search probably missed. Reformulate:
  try synonyms, the exact CLI flag, the verbatim error message.
- **If two reformulations both miss:** the corpus may not cover this
  question. See "Anti-hallucination" below.

### Step 3 — Fetch with the right `mode`

`fetch_doc(url_or_path, source, mode=...)` projections:

| `mode` | Use for |
|---|---|
| `markdown` (default) | Paper sections, short forum threads, short code files. Full body. |
| `outline` | **Scout first** on long docs pages, long forum threads, unfamiliar code files. Headings (docs/paper), defs/classes (code), per-post summary (forum). |
| `sections:<heading>` | Pull one named section from a docs page or paper section. Case-insensitive. |
| `post:<N>` | Forum only. One specific post's body, with a per-post URL. |
| `post:accepted` | Forum only. The accepted-answer post. Cheap and high-signal for solved threads. Errors if the thread is unsolved — use `outline` first if uncertain. |

Default to the cheaper projection (`outline`, `sections:…`,
`post:accepted`) when you can. Reach for `markdown` when you actually
need the full body.

### Step 4 — Cross-source when it helps

Some questions are best answered by combining sources:

- **"Why does X work this way and how do I use it?"** → paper for the
  why, docs for the how.
- **"I'm getting this error from HybridNew"** → forum for the
  diagnosis, docs for HybridNew context.
- **"What does --robustHesse actually do under the hood?"** → docs for
  the description, code for the implementation.

Don't query all four sources routinely — that wastes calls. Default
to one source per question. Parallel-querying two sources upfront is
fine when the question genuinely spans dimensions (like "how AND
why", or "reproduce the error AND explain the cause"). Fall through
to a third source only when the first two leave a real gap.

## Output format

- Start with the **direct answer** in 1–3 sentences.
- Cite the URLs returned by the tools, **verbatim**. Use inline
  Markdown links. Combine users will click through to verify.
- For forum citations, use the per-post URL when you fetched a
  specific post (`post:N` or `post:accepted`); the topic URL otherwise.
- If multiple sources contributed, list all relevant citations.
- Use code blocks for actual code or CLI invocations only — not for
  prose.

## Anti-hallucination

Combine is a domain where wrong answers can lead to wrong physics
results. Be conservative:

- If you can't find the answer in the corpus after two reformulations,
  **say so explicitly**: "I couldn't find this in the Combine docs,
  paper, code, or forum."
- Suggest the user post on
  [cms-talk](https://cms-talk.web.cern.ch/c/physics/cat/cat-stats/279).
- Do **not** answer from prior knowledge without flagging it. If you
  add context beyond what the corpus returned (e.g. comparing a
  fetched method to a related one, or noting a well-known
  consequence), say so explicitly: "not covered in the fetched
  section, but…" or "this is well-known but not in the corpus".
- Never paraphrase a forum reply as if it were the canonical doc.
  Cite the thread and let the user judge.

## Worked example (retrieval)

> *User:* "I'm running AsymptoticLimits and getting 'cannot compute
> the expected limit'. What's wrong?"

Reasoning: this is an error message — start with the forum.

1. `search_docs(query="cannot compute expected limit AsymptoticLimits", source="combine-forum")`
   → look at the top hit(s).
2. If a solved thread surfaces: `fetch_doc(url_or_path="<topic_id>", source="combine-forum", mode="post:accepted")`.
3. If the thread is unsolved:
   `fetch_doc(url_or_path="<topic_id>", source="combine-forum", mode="outline")`,
   then fetch the most-relevant reply with `post:<N>`.
4. Optional context:
   `search_docs(query="AsymptoticLimits", source="combine-docs")` →
   fetch the section that explains what the expected limit
   calculation does.
5. Answer: explain the cause (cite the cms-talk URL), point at the
   fix (cite the post URL), and link the docs page for canonical
   reference.

---

# Part 2 — Running Combine (execution)

Before running anything, decide **how** to execute — and the rule is
simple: **if `combine` is on the user's `PATH`, run it in the shell;
otherwise use the `combine-run-remote` server.** Check with
`command -v combine` first, every time (see "Choosing how to execute").
Getting this wrong is the most common failure — do not default to the
`run_combine` tool just because it appears in your tool list.

The two execution paths, in order of preference:

- **Shell** (Combine on `PATH`) — the default. The user sourced their
  environment (e.g. `cmsenv`) before launching the agent, so commands
  run in their real working directory against their real files, by path.
- **`combine-run-remote`** server — only when Combine is *not* on
  `PATH`. Runs each command in an isolated, throwaway workspace, so you
  pass every input file explicitly.

In both cases, run **one** Combine-family command at a time (`combine`,
`text2workspace.py`, `combineCards.py`, `combineTool.py`) — no shell
pipes, redirects, or chained commands.

## When to run vs. when to just explain

**Diagnose from the command first.** For a *reported failing command*,
check whether the cause is apparent from the command itself — an option
set outside its allowed range, mutually inconsistent flags, a typo, a
wrong file or mode. If it is, explain the cause and the fix **directly**
(cross-checking the `combine` sources for the exact rule); running it is
then only optional confirmation. Reproduce by running when the cause is
*not* apparent from inspection, or when you need the real error text.

**Run** when:
- The user reports a command that errors or gives an unexpected result,
  the cause is **not** apparent from the command itself, and you can
  reproduce it — running is then the fastest path to a real diagnosis.
- The user explicitly asks you to run, try, or check something.
- You've proposed a fix and want to confirm it works.

**Do not run — explain from the corpus instead** when:
- **The failure is explainable from the command/options alone** — e.g. a
  parameter set outside its declared range, mutually inconsistent flags,
  or an obvious typo. Explain the cause and the fix directly; running is
  optional confirmation, not the first move.
- **No execution is available** — Combine isn't on `PATH` and the
  `combine-run-remote` server isn't registered. Never fabricate
  execution results. Say execution isn't available, then reason about the
  command from the `combine` sources (what the flag does, what the error
  usually means).
- The task is batch submission (e.g. `combineTool.py --job-mode condor`)
  — don't submit jobs on the user's behalf; explain instead.
- The inputs are large (see routing below) and only the size-capped
  `combine-run-remote` server is available (no Combine on `PATH`).

## Choosing how to execute

**Always do this check FIRST, before every execution — do not skip it
and do not assume the answer.** Run `command -v combine` in the shell
(your bash/shell tool). Its result decides everything below:

- **`combine` IS on `PATH`** → you **must** run in the **shell**. Do
  **not** call the `combine-run-remote` tool in this case, even though it
  is registered and visible in your tool list. The user has Combine
  installed and their files on disk; running in the shell uses their
  real files in place, by path — which is what they expect. Reaching for
  `run_combine` here is a bug: it runs in an isolated throwaway workspace
  that cannot see the user's files, so it fails to find them.
- **`combine` is NOT on `PATH`** → use the `combine-run-remote` server.

The mere existence of the `combine-run-remote` tool is **not** a reason
to use it. It is a *fallback* for when the user has no local Combine.
The shell is the default whenever Combine is on `PATH`. When
`command -v combine` prints a path, `combine-run-remote` is **off the
table for this session** — treat it as if it did not exist.

> **Concrete example of the mistake to avoid.**
> User: "run `combine -M AsymptoticLimits -d data/tutorials/counting/realistic-counting-experiment.txt`".
>
> WRONG — what not to do:
> `command -v combine` → `/…/CMSSW_14_1_0_pre4/…/bin/combine` (Combine
> IS on PATH), yet the model calls
> `combine-run-remote_run_combine(command="combine -M AsymptoticLimits -d data/tutorials/…")`.
> It fails: the `run_combine` workspace is isolated and empty, so the
> relative path `data/tutorials/…` isn't there. The model reached for
> the tool named "run combine" instead of the shell — the exact bug.
>
> RIGHT — what to do:
> `command -v combine` prints a path → run the command with your **shell
> tool**, in the user's working directory, so the relative path resolves
> against their real files:
> `combine -M AsymptoticLimits -d data/tutorials/counting/realistic-counting-experiment.txt`
> Then read stdout for the expected limit. No `run_combine` call at all.

The two ways of "running combine" are genuinely different tools: the
**shell tool** runs the real executable in the user's directory against
their files; the **`combine-run-remote` MCP tool** copies files into a
throwaway sandbox on a remote server. When the user has Combine locally,
only the shell does what they mean.

Order of preference:

1. **Shell** — Combine on `PATH`. Runs against the user's real files in
   their working directory, no size caps, outputs persist. The best
   default. **(See "Running in the shell" below.)**
2. **`combine-run-remote`** — the shared CERN PaaS service, used only
   when Combine is *not* on `PATH`. Isolated workspace, size-capped
   inputs, shorter timeouts; you pass input files explicitly (see
   "Running via the `combine-run-remote` server").

If neither is available, don't run — explain from the corpus.

Note — **large inputs** (datacards + shape files more than a few MB):
only the shell handles these; `combine-run-remote` will reject them, so
if Combine isn't on `PATH` and the inputs are large, explain rather than
run.

## Running in the shell

When Combine is on `PATH`:

- Run one command, e.g.
  `combine -M AsymptoticLimits -d datacard.txt -m 125`. It executes in
  the user's current working directory, so the datacard and any shape
  files it references are expected to be on disk there already — don't
  invent file contents.
- `combine`, `text2workspace.py`, and `combineCards.py` run without a
  prompt (they're allow-listed). `combineTool.py` prompts, because it
  can submit batch jobs — and you should **not** submit batch jobs
  (condor/crab) on the user's behalf; explain instead.
- Read the exit code, stdout, and stderr directly. Output files (e.g.
  `higgsCombine*.root`) land in the working directory.

## Running via the `combine-run-remote` server

- `command`: the full command line, e.g.
  `"combine -M AsymptoticLimits -d datacard.txt -m 125"`. One command;
  no shell pipes or redirects.
- `files`: text inputs as `{filename: content}` — the datacard, any
  text configs. Reference them in the command by the same filename.
- `files_b64`: binary inputs as `{filename: base64}` — ROOT shape
  files a datacard references.
- `timeout_s`: optional; the server clamps it to its ceiling.

The server runs the command in an isolated, throwaway workspace on a
remote machine: **only the files you pass exist there.** It does not see
the user's filesystem. So a path argument to a file you did not pass
will fail with "file not found".

**When the user gives you a path to a datacard, read it yourself — do
not ask them to paste its contents.** You have a file-reading tool; use
it:

1. Read the datacard at the given path.
2. Scan it for `shapes` lines — each names a ROOT file the datacard
   references. Read those too (as binary → base64).
3. Call `run_combine` passing the datacard via `files` and every
   referenced ROOT file via `files_b64`, and reference them in the
   `command` by their **base name** (e.g. `-d datacard.txt`, not the
   original absolute path — the workspace is flat).

Only ask the user for a file if you genuinely cannot read it (it does
not exist, or is outside a path you can access). Pass every file the
command needs: a shape-based datacard that references `shapes.root`
won't run unless you also pass that file via `files_b64`.

## Interpreting the result

In the shell, read the exit code, stdout, and stderr directly — a
non-zero exit is the debugging case (cross-check `stderr` against the
`combine` sources). The `combine-run-remote` server instead returns
JSON that distinguishes three outcomes:

- **`error` is set** (and `returncode` is null): a *setup* failure —
  disallowed executable, oversized input, unsafe filename, or Combine
  not found on PATH. This is not a physics result. If it says "not
  found on PATH", the server has no Combine environment — tell the
  user; don't guess at numbers.
- **`returncode` is non-zero**: Combine ran but failed. Read `stderr`.
  This is the interesting debugging case — cross-check the error
  against the `combine` sources (`combine-forum` for the message text,
  `combine-code` for what the code checks) to explain and fix.
- **`returncode` is 0**: success. Report the results from `stdout` and
  note the `artifacts` produced (e.g. `higgsCombine*.root`).

Watch `timed_out` (bump `timeout_s`, or if the job is too heavy for the
size-capped remote server, note that it needs a local Combine) and the
`*_truncated` flags (output was tailed).

## The reproduce → diagnose → fix loop

1. **Reproduce:** run the user's command — in the shell if Combine is on
   `PATH`, otherwise via `run_combine` with the files passed in.
2. **Read the outcome:** setup error vs. non-zero exit vs. success.
3. **Diagnose (if it failed):** cross-check `stderr` against the
   `combine` retrieval sources — forum for "has anyone hit this",
   code/docs for what the failing option requires.
4. **Fix:** propose a corrected command, explaining what changed.
5. **Confirm:** re-run the corrected command; report the new result
   with citations for the reasoning.

## Worked example (execution + diagnosis)

> *User:* "This errors: `combine -M FitDiagnostics -d card.txt
> --robustHesse 1`, and here's my datacard."

1. Run `combine -M FitDiagnostics -d card.txt --robustHesse 1` — in the
   shell if Combine is on `PATH` (the card is already on disk), otherwise
   `run_combine(command=..., files={"card.txt": <their card>})`.
2. Inspect the result: if Combine wasn't found, execution isn't
   available — say so, don't guess numbers. If it exited non-zero, read
   `stderr`.
3. Take the key line from `stderr` and
   `search_docs(query="<that error text>", source="combine-forum")`;
   if needed, check `combine-code` for what the failing step requires.
4. Propose the fix (cite the forum thread / docs), then re-run the
   corrected command to confirm.
5. Report: what failed, why (with citations), and the confirmed
   working command.
