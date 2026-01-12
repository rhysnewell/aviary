When adding new binning tools:
1. When running a specific test, use `pixi run run-a-test -- test_name`, where `test_name` is the name of the test function to run (it is specified as the argument for pytest's `-k`). This allows for quick iteration on a specific test without running the full test suite.
2. To test for changes more broadly (but still relatively quickly), use `pixi run run-quick-tests`
3. Add at least one test to test/test_integration.py which actually runs the tool being tested.
4. Add at least one test to test/test_recovery.py which makes sure that the intermediate config file is accurately generated.
5. Update the documentation in docs/binning_tools.md to include usage instructions and examples for the new tool. 
6. Update docs/citations.md to include any relevant citations for the new tool.
