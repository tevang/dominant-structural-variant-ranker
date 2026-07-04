.PHONY: ci-log

ci-log:
	@test -n "$(RUN)" || (echo "Usage: make ci-log RUN=<actions-run-url-or-id> [REPO=owner/repo]" && exit 1)
	python scripts/inspect_ci_run.py $(RUN) $(if $(REPO),--repo $(REPO),)
