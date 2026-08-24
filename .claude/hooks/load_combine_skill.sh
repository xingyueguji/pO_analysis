#!/bin/bash
# SessionStart hook: load the project-scoped `combine` skill into context.
#
# Why a hook: the skill is discoverable in every session anyway (it lives in
# .claude/skills/), but discoverable only means Claude sees its one-line
# description and may or may not invoke it. This injects the actual content, so
# Combine knowledge is present from turn one.
#
# MODE=full     inject the whole SKILL.md (default)
# MODE=pointer  inject only a short directive to invoke the skill on demand
#               (cheap; use this if the full text costs too much context)
MODE="${COMBINE_SKILL_MODE:-full}"

set -u

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
skill="$here/../skills/combine/SKILL.md"

[ -r "$skill" ] || exit 0   # skill removed/renamed -> stay silent, don't break startup

MODE="$MODE" SKILL="$skill" python3 <<'PY'
import json, os

mode  = os.environ["MODE"]
path  = os.environ["SKILL"]

# Status of the retrieval backend this skill was written against. The `combine`
# MCP server (combine-mcp-git-combine-mcp.app.cern.ch) was removed 2026-08-20
# after it stopped responding: DNS + TLS succeed, then the server closes the
# connection with no HTTP reply on every path. Until it is restored, Part 1's
# search_docs/fetch_doc tools DO NOT EXIST -- do not attempt to call them.
status = (
    "STATUS (2026-08-20): the `combine` MCP server is NOT configured -- the CERN "
    "endpoint it pointed at is down. The skill's Part 1 (retrieval) tools "
    "`search_docs`/`fetch_doc` are therefore unavailable; do not call them. "
    "Answer Combine questions from this repo's own docs (CLAUDE.md \"Downstream "
    "fit\", the fork's test/README_pO_fits.md) and from the shell-execution half "
    "of Part 2. Re-add the server to .mcp.json if it comes back."
)

if mode == "pointer":
    body = (
        "A project-scoped `combine` skill covers CMS Combine (datacards, physics "
        "models, running modes, error messages, interpreting output). Invoke it "
        f"with the Skill tool (skill: \"combine\") for any Combine work; its full "
        f"text is at {path}.\n\n{status}"
    )
else:
    with open(path, encoding="utf-8") as fh:
        body = fh.read()
    body = (
        "The project-scoped `combine` skill is loaded below (auto-injected at "
        f"session start). Source: {path}\n\n{status}\n\n---\n\n{body}"
    )

print(json.dumps({
    "suppressOutput": True,
    "hookSpecificOutput": {
        "hookEventName": "SessionStart",
        "additionalContext": body,
    },
}))
PY
