hw03_worker = hw03();

hw_assert(abs(hw03_worker.p1([-1, 1;1, 1], 0) - 1) < 1e-8);
hw_assert(abs(hw03_worker.p1([-1, 1;1, 2], 0) - 1.5)<1e-8);
hw_assert(abs(hw03_worker.p1([-1, 1;1, 2;3, 3], 0) - 1.5)<1e-8);
hw_assert(abs(hw03_worker.p1([-1, 1;1, 2;3, 7], 0) - 1)<1e-8);
hw_assert(abs(hw03_worker.p1([-1, 1;1, 2;3, 7; 4, 2], 0) + 0.8)<1e-8);
hw_assert(abs(hw03_worker.p1([-1, 1;1, 2;3, 7; 4, 2; 5, 9], 0) + 6.25)<1e-8);


% hw03_worker.p2();
% hw03_worker.p3();

function hw_assert(X)
    if X; fprintf('\t PASS\n'); else; fprintf('\t FAIL\n'); end
end