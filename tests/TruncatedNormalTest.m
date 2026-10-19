classdef (SharedTestFixtures = {ProteusTestEnvironment}) ...
        TruncatedNormalTest < matlab.unittest.TestCase
%TRUNCATEDNORMALTEST random_truncated_normal (base MATLAB) must reproduce
%the Statistics and Machine Learning Toolbox draw used by the original
%define_medium.

    properties (TestParameter)
        seed   = {0, 3, 42};
        cutoff = {1, 2, 2.5, 3};
        sz     = {[7 11 13], [64 64 64]};
    end

    methods (Test, TestTags = {'Fast'})
        function matchesToolboxBitForBit(testCase, seed, cutoff, sz)
            testCase.assumeTrue(~isempty(ver('stats')) && ...
                license('test','Statistics_Toolbox'), ...
                'Statistics and Machine Learning Toolbox not available.');
            testCase.addTeardown(@rng, rng);

            pd = truncate(makedist('Normal'), -cutoff, cutoff);
            rng(seed,'twister'); expected = random(pd, sz);
            rng(seed,'twister'); actual   = random_truncated_normal(cutoff, sz);

            testCase.verifyEqual(actual, expected);
            % Same number of draws from the global stream:
            rng(seed,'twister'); random(pd, sz); next1 = rand;
            rng(seed,'twister'); random_truncated_normal(cutoff, sz); next2 = rand;
            testCase.verifyEqual(next2, next1);
        end

        function infiniteCutoffIsUntruncated(testCase)
            testCase.addTeardown(@rng, rng);
            rng(5,'twister'); expected = randn([4 5 6]);
            rng(5,'twister'); actual = random_truncated_normal(Inf, [4 5 6]);
            testCase.verifyEqual(actual, expected);
        end

        function invalidCutoffErrors(testCase)
            for c = {0, -1, NaN, [1 2]}
                testCase.verifyError(@() random_truncated_normal(c{1}, [2 2]), ...
                    'random_truncated_normal:invalidCutoff');
            end
        end

        function staysWithinCutoff(testCase, cutoff)
            testCase.addTeardown(@rng, rng);
            rng(1,'twister');
            x = random_truncated_normal(cutoff, [100 100 10]);
            testCase.verifyLessThanOrEqual(max(abs(x(:))), cutoff);
            testCase.verifyEqual(mean(x(:)), 0, 'AbsTol', 0.02);
        end
    end
end
