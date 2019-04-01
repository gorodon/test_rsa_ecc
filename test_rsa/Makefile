
BIN_NAME = test_rsa

OBJS= test.o \
      rsaimpl.o \
      ../Vlong.o \
      ../Vlong_Gp.o

ifeq ("$(RELEASEDIR)","")
  RELEASEDIR=../nix_release
endif

CFLAGS +=-std=c99
CFLAGS += -Ofast
CFLAGS += -Wall
CFLAGS += -I. -I../

#hidden rules
#.c.o:
#       $(CC) $(CFLAGS) -c $<

.PHONY: all
all: $(BIN_NAME)

$(BIN_NAME): $(OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LDFLAGS2)
	mkdir -p $(RELEASEDIR)
	cp -f $(BIN_NAME) $(RELEASEDIR)
	@echo "  HOST    Done compiling $(BIN_NAME)"

.cpp.o:
	$(CXX) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	@rm -f $(OBJS) $(BIN_NAME)
