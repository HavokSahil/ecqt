CC			:= gcc
CFLAGS 		:= -g -O0 -Wall -lm
INC			:= .

CFLAGS 		+= $(INC:%=-I%)

.PHONY: default
default: all

pffft.o: pffft/pffft.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: tests/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test_sm%: test_sm%.o
	$(CC) $(CFLAGS) $< -o $@

test_cqt%: test_cqt%.o pffft.o
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: test_sm
test_sm: 	test_sm_basic \
			test_sm_complex \
			test_sm_memory \
			test_sm_performance
	@echo "Performing all Sparse Matrix tests..."
	@./test_sm_basic
	@./test_sm_complex
	@./test_sm_memory
	@./test_sm_performance
	@echo "All Sparse Matrix tests finished."

.PHONY: test_cqt
test_cqt: 	test_cqt_basic \
			test_cqt_memory \
			test_cqt_performance \
			test_cqt_sparsekernel \
			test_cqt_transform
	@echo "Performing all CQT tests..."
	@./test_cqt_basic
	@./test_cqt_memory
	@./test_cqt_performance
	@./test_cqt_sparsekernel
	@./test_cqt_transform
	@echo "All CQT tests finished."

.PHONY: all
all: test_sm test_cqt

.PHONY: clean
clean:
	@rm -f \
		tests/*.o \
		*.o \
		test_*
