function tests = maskOrderingTest
addpath("../src")
tests = functiontests(localfunctions);
end

function testSimpleOrdering(testCase)
    v = [1 0 0 1 0 1 1 0 0 1];
    P = getMaskOrderingMatrix(v);
    verifyEqual(testCase,P*v',[1 1 1 1 1 0 0 0 0 0]');
end

function testTimeTaken(testCase)
    tic
    v = [zeros(200000,1) ones(200000,1)];
    P = getMaskOrderingMatrix(v);
    time = toc;
    fprintf("Time to order mask matrix: %.2f",time);
    verifyLessThan(testCase,time,10);
end