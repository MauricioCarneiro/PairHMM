CODE_DIR = src

.PHONY: project_code

project_code:
	$(MAKE) -C $(CODE_DIR)

clean:
	$(MAKE) -C $(CODE_DIR) clean

check:
	valgrind --leak-check=yes ./pairhmm test_data/tiny.in > /dev/null
