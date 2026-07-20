#!/usr/bin/env bash
#
# Read-only initial triage for Documenter.jl deployment credentials in a
# GitHub organization. Secret values are never available to this script; it
# only reads secret names/timestamps and public deploy-key metadata.

set -uo pipefail

org="${GITHUB_ORG:-QuantumSavory}"
format="tsv"
include_archived=true
include_forks=false

usage() {
    cat <<'EOF'
Usage: scripts/audit_documenter_setup.sh [options]

Audit Julia packages in a GitHub organization for Documenter.jl deployment
metadata. Requires authenticated `gh` and `jq`. Repository admin access is
needed for definitive secret and deploy-key results.

Options:
  --org ORG             GitHub organization (default: QuantumSavory)
  --format FORMAT       Output format: tsv or jsonl (default: tsv)
  --active-only         Exclude archived repositories
  --include-forks       Include repositories forked into the organization
  -h, --help            Show this help

The script deliberately keeps two conclusions separate:

  ssh_setup
      Structural metadata only. METADATA_PRESENT means an effective Actions
      secret named DOCUMENTER_KEY and a plausible verified writable deploy key
      both exist. GitHub does not expose the secret value, so the script cannot
      prove that the private and public halves match.

  triage
      The likely operational class. TOKEN_FALLBACK_AVAILABLE means the classic
      SSH setup is absent but GITHUB_TOKEN has apparent write permission. That
      can deploy same-repository dev docs, while TagBot-created tags still need
      SSH to trigger downstream versioned-doc workflows.

Important triage values include:
  METADATA_PRESENT, TOKEN_FALLBACK_AVAILABLE, SECRET_WITHOUT_WRITE_KEY,
  DEPLOY_TARGET_MISMATCH, DUPLICATE_DOC_DEPLOY, BUILDKITE_DOCS,
  DEAD_DOC_JOB_NO_MAKE, DOCS_BUILD_ONLY, NO_DEPLOYDOCS,
  OUT_OF_REPO_DEPLOY, and METADATA_UNAVAILABLE.
EOF
}

while (($# > 0)); do
    case "$1" in
        --org)
            if (($# < 2)); then
                printf 'error: --org requires a value\n' >&2
                exit 2
            fi
            org="$2"
            shift 2
            ;;
        --format)
            if (($# < 2)); then
                printf 'error: --format requires a value\n' >&2
                exit 2
            fi
            format="$2"
            shift 2
            ;;
        --active-only)
            include_archived=false
            shift
            ;;
        --include-forks)
            include_forks=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            printf 'error: unknown argument: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ "$format" != "tsv" && "$format" != "jsonl" ]]; then
    printf 'error: --format must be tsv or jsonl\n' >&2
    exit 2
fi

for command_name in gh jq; do
    if ! command -v "$command_name" >/dev/null 2>&1; then
        printf 'error: required command not found: %s\n' "$command_name" >&2
        exit 2
    fi
done

if ! gh auth status >/dev/null 2>&1; then
    printf 'error: gh is not authenticated\n' >&2
    exit 2
fi

raw_file() {
    local repo="$1"
    local path="$2"
    local ref="$3"
    local encoded_ref
    encoded_ref="$(jq -nr --arg value "$ref" '$value | @uri')"
    gh api "repos/$repo/contents/$path?ref=$encoded_ref" \
        -H 'Accept: application/vnd.github.raw+json' 2>/dev/null
}

matches() {
    local pattern="$1"
    local text="$2"
    grep -Eqi "$pattern" <<<"$text"
}

without_comment_lines() {
    sed '/^[[:space:]]*#/d'
}

extract_quoted_assignment() {
    local assignment="$1"
    local text="$2"
    local regex="(^|[^[:alnum:]_])${assignment}[[:space:]]*=[[:space:]]*\"([^\"]+)\""
    if [[ "$text" =~ $regex ]]; then
        printf '%s' "${BASH_REMATCH[2]}"
        return 0
    fi
    return 1
}

normalize_github_repo() {
    local value="$1"
    value="${value#ssh://git@github.com/}"
    value="${value#git@github.com:}"
    value="${value#https://github.com/}"
    value="${value#http://github.com/}"
    value="${value#github.com/}"
    value="${value%/}"
    value="${value%.git}"
    printf '%s' "$value"
}

repo_fields='nameWithOwner,name,isArchived,isFork,isPrivate,defaultBranchRef,pushedAt,stargazerCount,url'
if ! repos_json="$(gh repo list "$org" --limit 1000 --json "$repo_fields")"; then
    printf 'error: failed to list repositories for %s\n' "$org" >&2
    exit 1
fi

if [[ "$format" == "tsv" ]]; then
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "repository" "archived" "stars" "pushed_at" "deploydocs" \
        "deploy_target" "target_match" "docs_ci" "docs_workflows" \
        "docs_auth" "token_permission" "tagbot_ssh" "secret_origin" \
        "deploy_key" "ssh_setup" "triage"
fi

package_count=0
credential_candidate_count=0

while IFS= read -r repo_json; do
    repo="$(jq -r '.nameWithOwner' <<<"$repo_json")"
    archived="$(jq -r '.isArchived' <<<"$repo_json")"
    fork="$(jq -r '.isFork' <<<"$repo_json")"
    branch="$(jq -r '.defaultBranchRef.name // ""' <<<"$repo_json")"
    stars="$(jq -r '.stargazerCount' <<<"$repo_json")"
    pushed_at="$(jq -r '.pushedAt // ""' <<<"$repo_json")"

    if [[ "$include_forks" != true && "$fork" == true ]]; then
        continue
    fi
    if [[ "$include_archived" != true && "$archived" == true ]]; then
        continue
    fi
    if [[ -z "$branch" ]]; then
        continue
    fi

    encoded_branch="$(jq -nr --arg value "$branch" '$value | @uri')"
    if ! tree_json="$(gh api "repos/$repo/git/trees/$encoded_branch?recursive=1" \
        2>/dev/null)"; then
        printf 'warning: could not inspect tree for %s\n' "$repo" >&2
        continue
    fi
    if [[ "$(jq -r '.truncated // false' <<<"$tree_json")" == true ]]; then
        printf 'warning: recursive tree was truncated for %s\n' "$repo" >&2
    fi

    root_project="$(jq -r \
        'any(.tree[]?; .type == "blob" and .path == "Project.toml")' \
        <<<"$tree_json")"
    root_module="$(jq -r \
        'any(.tree[]?; .type == "blob" and (.path | test("^src/[^/]+\\.jl$")))' \
        <<<"$tree_json")"
    if [[ "$root_project" != true || "$root_module" != true ]]; then
        continue
    fi

    if ! project_toml="$(raw_file "$repo" "Project.toml" "$branch")"; then
        printf 'warning: could not read Project.toml for %s\n' "$repo" >&2
        continue
    fi
    if ! matches '^[[:space:]]*name[[:space:]]*=' "$project_toml" ||
        ! matches '^[[:space:]]*uuid[[:space:]]*=' "$project_toml"; then
        continue
    fi
    ((package_count += 1))

    has_docs_make="$(jq -r \
        'any(.tree[]?; .type == "blob" and .path == "docs/make.jl")' \
        <<<"$tree_json")"
    docs_make=""
    if [[ "$has_docs_make" == true ]]; then
        docs_make="$(raw_file "$repo" "docs/make.jl" "$branch" || true)"
    fi
    docs_make_active="$(without_comment_lines <<<"$docs_make")"

    deploydocs=false
    if matches '(^|[^[:alnum:]_])deploydocs[[:space:]]*\(' \
        "$docs_make_active"; then
        deploydocs=true
    fi

    source_target="$(extract_quoted_assignment "repo" \
        "$docs_make_active" || true)"
    deploy_repo_target="$(extract_quoted_assignment "deploy_repo" \
        "$docs_make_active" || true)"
    out_of_repo=false
    if [[ -n "$deploy_repo_target" ]]; then
        deploy_target="$deploy_repo_target"
        out_of_repo=true
    else
        deploy_target="$source_target"
    fi

    target_match="not-applicable"
    normalized_target=""
    if [[ "$deploydocs" == true ]]; then
        if [[ -n "$deploy_target" ]]; then
            normalized_target="$(normalize_github_repo "$deploy_target")"
            if [[ "$out_of_repo" == true ]]; then
                target_match="external"
            elif [[ "${normalized_target,,}" == "${repo,,}" ]]; then
                target_match="true"
            else
                target_match="false"
            fi
        else
            target_match="unknown"
        fi
    fi

    workflows_api_available=false
    workflows_json='{"workflows":[]}'
    if workflows_json="$(gh api "repos/$repo/actions/workflows?per_page=100" \
        2>/dev/null)"; then
        workflows_api_available=true
    fi

    gha_docs=false
    gha_doc_key_ref=false
    gha_doc_token_ref=false
    gha_explicit_contents_write=false
    tagbot_key_ref=false
    gha_docs_count=0
    gha_doc_paths=()

    while IFS= read -r workflow_path; do
        [[ -z "$workflow_path" ]] && continue

        workflow_state="unknown"
        if [[ "$workflows_api_available" == true ]]; then
            workflow_state="$(jq -r --arg path "$workflow_path" \
                '[.workflows[] | select(.path == $path) | .state][0] // "unknown"' \
                <<<"$workflows_json")"
            if [[ "$workflow_state" != "active" ]]; then
                continue
            fi
        fi

        workflow="$(raw_file "$repo" "$workflow_path" "$branch" || true)"
        workflow_active="$(without_comment_lines <<<"$workflow")"
        is_docs_workflow=false
        if matches 'julia-actions/julia-docdeploy|docs/make\.jl' \
            "$workflow_active"; then
            gha_docs=true
            is_docs_workflow=true
            ((gha_docs_count += 1))
            gha_doc_paths+=("$workflow_path")
        fi
        if [[ "$is_docs_workflow" == true ]] &&
            matches 'DOCUMENTER_KEY' "$workflow_active"; then
            gha_doc_key_ref=true
        fi
        if [[ "$is_docs_workflow" == true ]] &&
            matches 'GITHUB_TOKEN' "$workflow_active"; then
            gha_doc_token_ref=true
        fi
        if [[ "$is_docs_workflow" == true ]] &&
            matches 'contents:[[:space:]]*write' "$workflow_active"; then
            gha_explicit_contents_write=true
        fi
        if matches 'JuliaRegistries/TagBot' "$workflow_active" &&
            matches 'DOCUMENTER_KEY' "$workflow_active"; then
            tagbot_key_ref=true
        fi
    done < <(
        jq -r '
            .tree[]?
            | select(.type == "blob")
            | .path
            | select(test("^\\.github/workflows/[^/]+\\.ya?ml$"))
        ' <<<"$tree_json"
    )

    buildkite_docs=false
    while IFS= read -r buildkite_path; do
        [[ -z "$buildkite_path" ]] && continue
        pipeline="$(raw_file "$repo" "$buildkite_path" "$branch" || true)"
        pipeline_active="$(without_comment_lines <<<"$pipeline")"
        if matches \
            'docs/make\.jl|Documenter\.Buildkite|julia-actions/julia-docdeploy' \
            "$pipeline_active"; then
            buildkite_docs=true
        fi
    done < <(
        jq -r '
            .tree[]?
            | select(.type == "blob")
            | .path
            | select(
                test("^\\.buildkite/.*\\.ya?ml$")
                or test("^buildkite\\.ya?ml$")
            )
        ' <<<"$tree_json"
    )

    if [[ "$gha_docs" == true && "$buildkite_docs" == true ]]; then
        docs_ci="github-actions+buildkite"
    elif [[ "$gha_docs" == true ]]; then
        docs_ci="github-actions"
    elif [[ "$buildkite_docs" == true ]]; then
        docs_ci="buildkite"
    else
        docs_ci="none"
    fi

    if [[ "$gha_doc_key_ref" == true && "$gha_doc_token_ref" == true ]]; then
        docs_auth="documenter-key+github-token"
    elif [[ "$gha_doc_key_ref" == true ]]; then
        docs_auth="documenter-key"
    elif [[ "$gha_doc_token_ref" == true ]]; then
        docs_auth="github-token"
    else
        docs_auth="none"
    fi

    token_permission="not-checked"
    if [[ "$gha_doc_token_ref" == true ]]; then
        if [[ "$gha_explicit_contents_write" == true ]]; then
            token_permission="write-explicit"
        elif token_permission_json="$(gh api \
            "repos/$repo/actions/permissions/workflow" 2>/dev/null)"; then
            default_token_permission="$(jq -r \
                '.default_workflow_permissions // "unknown"' \
                <<<"$token_permission_json")"
            token_permission="${default_token_permission}-default"
        else
            token_permission="unknown"
        fi
    fi

    secret_origin="not-checked"
    secret_updated_at=""
    key_state="not-checked"
    key_display="not-checked"
    ssh_setup="NOT_APPLICABLE"

    if [[ "$gha_docs" == true && "$deploydocs" == true &&
        "$out_of_repo" != true ]]; then
        ((credential_candidate_count += 1))

        repo_secret_state="unknown"
        org_secret_state="unknown"
        repo_secret_updated=""
        org_secret_updated=""

        if repo_secrets_json="$(gh api "repos/$repo/actions/secrets" \
            2>/dev/null)"; then
            if jq -e 'any(.secrets[]?; .name == "DOCUMENTER_KEY")' \
                <<<"$repo_secrets_json" >/dev/null; then
                repo_secret_state="present"
                repo_secret_updated="$(jq -r \
                    '.secrets[]
                     | select(.name == "DOCUMENTER_KEY")
                     | .updated_at' <<<"$repo_secrets_json")"
            else
                repo_secret_state="absent"
            fi
        fi

        if org_secrets_json="$(gh api \
            "repos/$repo/actions/organization-secrets" 2>/dev/null)"; then
            if jq -e 'any(.secrets[]?; .name == "DOCUMENTER_KEY")' \
                <<<"$org_secrets_json" >/dev/null; then
                org_secret_state="present"
                org_secret_updated="$(jq -r \
                    '.secrets[]
                     | select(.name == "DOCUMENTER_KEY")
                     | .updated_at' <<<"$org_secrets_json")"
            else
                org_secret_state="absent"
            fi
        fi

        if [[ "$repo_secret_state" == "present" ]]; then
            secret_origin="repository"
            secret_updated_at="$repo_secret_updated"
        elif [[ "$org_secret_state" == "present" ]]; then
            secret_origin="organization"
            secret_updated_at="$org_secret_updated"
        elif [[ "$repo_secret_state" == "unknown" ||
            "$org_secret_state" == "unknown" ]]; then
            secret_origin="unknown"
        else
            secret_origin="absent"
        fi

        doc_write_key_count=0
        doc_read_key_count=0
        any_write_key_count=0
        if keys_json="$(gh api "repos/$repo/keys" 2>/dev/null)"; then
            doc_write_key_count="$(jq '
                [
                    .[]
                    | select(.read_only == false and .verified == true)
                    | select(
                        ((.title // "") | test("documenter|documentation|docs"; "i"))
                        or ((.key // "") | test("\\sDocumenter$"; "i"))
                    )
                ]
                | length
            ' <<<"$keys_json")"
            doc_read_key_count="$(jq '
                [
                    .[]
                    | select(.read_only == true)
                    | select(
                        ((.title // "") | test("documenter|documentation|docs"; "i"))
                        or ((.key // "") | test("\\sDocumenter$"; "i"))
                    )
                ]
                | length
            ' <<<"$keys_json")"
            any_write_key_count="$(jq \
                '[.[] | select(.read_only == false and .verified == true)]
                 | length' <<<"$keys_json")"
            key_display="$(jq -r '
                [
                    .[]
                    | "\(.title // "(untitled)")[" +
                      (if .read_only then "read" else "write" end) + "]"
                ]
                | if length == 0 then "absent" else join(",") end
            ' <<<"$keys_json")"

            if ((doc_write_key_count > 0)); then
                key_state="documenter-write"
            elif ((doc_read_key_count > 0)); then
                key_state="documenter-read-only"
            elif ((any_write_key_count > 0)); then
                key_state="write-candidate"
            else
                key_state="absent"
            fi
        else
            key_state="unknown"
            key_display="unknown"
        fi

        if [[ "$secret_origin" == "unknown" || "$key_state" == "unknown" ]]; then
            ssh_setup="METADATA_UNAVAILABLE"
        elif [[ "$secret_origin" == "absent" && "$key_state" == "absent" ]]; then
            ssh_setup="MISSING_BOTH"
        elif [[ "$secret_origin" == "absent" ]]; then
            ssh_setup="MISSING_SECRET"
        elif [[ "$key_state" == "documenter-read-only" ]]; then
            ssh_setup="DEPLOY_KEY_READ_ONLY"
        elif [[ "$key_state" == "write-candidate" ]]; then
            ssh_setup="WRITE_KEY_AMBIGUOUS"
        elif [[ "$key_state" == "absent" ]]; then
            ssh_setup="MISSING_DEPLOY_KEY"
        else
            ssh_setup="METADATA_PRESENT"
        fi
    fi

    if [[ "$docs_ci" == "github-actions+buildkite" ]]; then
        triage="MIXED_DOCS_CI_MANUAL"
    elif [[ "$out_of_repo" == true ]]; then
        triage="OUT_OF_REPO_DEPLOY"
    elif [[ "$deploydocs" == true && "$target_match" == "false" ]]; then
        triage="DEPLOY_TARGET_MISMATCH"
    elif [[ "$docs_ci" == "buildkite" ]]; then
        triage="BUILDKITE_DOCS"
    elif [[ "$gha_docs" == true && "$has_docs_make" != true ]]; then
        triage="DEAD_DOC_JOB_NO_MAKE"
    elif [[ "$has_docs_make" == true && "$deploydocs" != true ]]; then
        triage="DOCS_BUILD_ONLY"
    elif [[ "$deploydocs" != true ]]; then
        triage="NO_DEPLOYDOCS"
    elif [[ "$gha_docs" != true ]]; then
        triage="DEPLOYDOCS_NO_SUPPORTED_CI"
    elif ((gha_docs_count > 1)); then
        triage="DUPLICATE_DOC_DEPLOY"
    elif [[ "$gha_doc_key_ref" != true && "$gha_doc_token_ref" == true ]]; then
        triage="GITHUB_TOKEN_ONLY"
    elif [[ "$gha_doc_key_ref" != true ]]; then
        triage="DOCUMENTER_KEY_NOT_EXPOSED"
    elif [[ "$ssh_setup" == "METADATA_PRESENT" ]]; then
        triage="METADATA_PRESENT"
    elif [[ "$ssh_setup" == "MISSING_DEPLOY_KEY" ||
        "$ssh_setup" == "DEPLOY_KEY_READ_ONLY" ]]; then
        triage="SECRET_WITHOUT_WRITE_KEY"
    elif [[ "$ssh_setup" == "METADATA_UNAVAILABLE" ]]; then
        triage="METADATA_UNAVAILABLE"
    elif [[ "$ssh_setup" == "WRITE_KEY_AMBIGUOUS" ]]; then
        triage="WRITE_KEY_MANUAL_REVIEW"
    elif [[ "$secret_origin" == "absent" &&
        "$gha_doc_token_ref" == true &&
        "$token_permission" == write-* ]]; then
        triage="TOKEN_FALLBACK_AVAILABLE"
    else
        triage="SSH_SETUP_INCOMPLETE"
    fi

    workflow_display="$(IFS=,; printf '%s' "${gha_doc_paths[*]}")"
    if [[ "$secret_origin" == "repository" ||
        "$secret_origin" == "organization" ]]; then
        secret_display="${secret_origin}@${secret_updated_at}"
    else
        secret_display="$secret_origin"
    fi

    if [[ "$format" == "tsv" ]]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$repo" "$archived" "$stars" "$pushed_at" "$deploydocs" \
            "${normalized_target:-}" "$target_match" "$docs_ci" \
            "$workflow_display" "$docs_auth" "$token_permission" \
            "$tagbot_key_ref" "$secret_display" "$key_display" \
            "$ssh_setup" "$triage"
    else
        workflow_paths_json="$(printf '%s\n' "${gha_doc_paths[@]}" |
            jq -Rsc 'split("\n") | map(select(length > 0))')"
        jq -nc \
            --arg repository "$repo" \
            --argjson archived "$archived" \
            --argjson stars "$stars" \
            --arg pushed_at "$pushed_at" \
            --argjson deploydocs "$deploydocs" \
            --arg deploy_target "$normalized_target" \
            --arg target_match "$target_match" \
            --arg docs_ci "$docs_ci" \
            --argjson docs_workflows "$workflow_paths_json" \
            --arg docs_auth "$docs_auth" \
            --arg token_permission "$token_permission" \
            --argjson tagbot_ssh "$tagbot_key_ref" \
            --arg secret_origin "$secret_origin" \
            --arg secret_updated_at "$secret_updated_at" \
            --arg deploy_key "$key_state" \
            --arg deploy_key_display "$key_display" \
            --arg ssh_setup "$ssh_setup" \
            --arg triage "$triage" \
            '{
                repository: $repository,
                archived: $archived,
                stars: $stars,
                pushed_at: $pushed_at,
                deploydocs: $deploydocs,
                deploy_target: $deploy_target,
                target_match: $target_match,
                docs_ci: $docs_ci,
                docs_workflows: $docs_workflows,
                docs_auth: $docs_auth,
                token_permission: $token_permission,
                tagbot_ssh: $tagbot_ssh,
                secret_origin: $secret_origin,
                secret_updated_at: $secret_updated_at,
                deploy_key: $deploy_key,
                deploy_key_display: $deploy_key_display,
                ssh_setup: $ssh_setup,
                triage: $triage
            }'
    fi
done < <(jq -c 'sort_by(.nameWithOwner)[]' <<<"$repos_json")

printf 'Audited %d Julia package repositories; checked SSH metadata for %d GitHub Actions deployers.\n' \
    "$package_count" "$credential_candidate_count" >&2
