EXEC = wt
ARGS = -nox -iprint 50 

	DIST = ../../src/

ifdef ComSpec
    RM=del /F /Q
else
    RM=rm -rf
endif

all: both

$(EXEC): $(DIST)$(EXEC).tpl
	ln -sf $(DIST)$(EXEC) $@

$(DIST)$(EXEC).tpl:
	$(MAKE) --directory=../../src

run: $(EXEC)
	./$(EXEC) $(ARGS)

plot:
	R CMD BATCH plot.R

both: $(EXEC)
	./$(EXEC) $(ARGS)
	R CMD BATCH plot.R

sim: $(EXEC)
	./$(EXEC) -ainp wt.par -sim 123

clean:
	@$(RM) -f $(EXEC)
	@$(RM) -f admodel.*
	@$(RM) -f $(EXEC)  $(EXEC).[brces]* $(EXEC).*[0123456789] *.rpt *.log variance

clobber: clean
	@$(RM) -f *.par *_in.* *.rep
