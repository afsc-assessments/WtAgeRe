
ifdef ComSpec
    RM=del /F /Q
    COPY=copy
else
    RM=rm -rf
    COPY=cp
endif

# Keep package installation portable: build ADMB binary only when available.
all:
ifdef ComSpec
	@where admb >NUL 2>NUL && (echo Building wt with ADMB... && admb wt.tpl || echo ADMB build failed; skipping optional wt build.) || (echo ADMB not found; skipping optional wt build.)
else
	@if command -v admb >/dev/null 2>&1; then \
		echo "Building wt with ADMB..."; \
		admb wt.tpl || echo "ADMB build failed; skipping optional wt build."; \
	else \
		echo "ADMB not found; skipping optional wt build."; \
	fi
endif

wt: wt.tpl 
ifdef ComSpec
	@where admb >NUL 2>NUL || (echo ADMB is required to build wt. && exit 1)
	@admb wt.tpl
else
	@command -v admb >/dev/null 2>&1 || { echo "ADMB is required to build wt." >&2; exit 1; }
	@admb wt.tpl
endif


clean:
	@$(RM) build
	@$(RM) wt.cpp
	@$(RM) wt.htp
	@$(RM) wt.obj
