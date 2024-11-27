// #include "interactionRewards.h"

// REWARD::REWARD() {

// }

// REWARD::REWARD(int timeSteps, double rewardStart, double rewardEnd, double weight) {
//     _ts = timeSteps;
//     _rs = rewardStart;
//     _re = rewardEnd;
//     _fixedLimit = rewardEnd;
//     _w = weight;
// }

// double REWARD::interactionReward(int step) {
//     //std::vector<int> rewards(steps, 0);
//     double rewards;
//     for (int i = _rs; i <= _fixedLimit; i++) {
//         double rewardCoefficient = static_cast<double>(i - _rs) / (_re - _rs);
//         rewards = rewardCoefficient * _w;
//     }

//     return rewards;
// }