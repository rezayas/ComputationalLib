SUBDIRS = lib test
all: $(SUBDIRS)
.PHONY: $(SUBDIRS) all clean
$(SUBDIRS):
	$(MAKE) -C $@
test: lib
clean:
	$(MAKE) clean -C lib
	$(MAKE) clean -C test
