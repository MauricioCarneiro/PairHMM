CODE_DIR = src

.PHONY: project_code

project_code:
	$(MAKE) -C $(CODE_DIR)

clean:
	$(MAKE) -C $(CODE_DIR) clean

test10s: project_code
	@./pairhmm < test_data/10s.in | paste - test_data/10s.out  | awk 'BEGIN {m = 0; n = NR} {a = $$2-$$1; a = a < 0? -a: a; if (a > m) {m = a; n = NR} } END { printf("max error %g at line %d\n", m, n)} '

test1m: project_code
	@./pairhmm < test_data/1m.in | paste - test_data/1m.out  | awk 'BEGIN {m = 0; n = NR} {a = $$2-$$1; a = a < 0? -a: a; if (a > m) {m = a; n = NR} } END { printf("max error %g at line %d\n", m, n)} '

testmedium: project_code
	@xzcat test_data/medium.in.xz | pairhmm | paste - test_data/medium.out  | awk 'BEGIN {m = 0; n = NR} {a = $$2-$$1; a = a < 0? -a: a; if (a > m) {m = a; n = NR} } END { printf("max error %g at line %d\n", m, n)} '

check: project_code
	valgrind --leak-check=yes ./pairhmm test_data/tiny.in > /dev/null
