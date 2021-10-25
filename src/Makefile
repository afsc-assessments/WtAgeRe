
ifdef ComSpec
    RM=del /F /Q
    COPY=copy
else
    RM=rm -rf
    COPY=cp
endif

all: wt

wt: wt.tpl 
	@admb wt.tpl 


clean:
	@$(RM) build
	@$(RM) wt.cpp
	@$(RM) wt.htp
	@$(RM) wt.obj
