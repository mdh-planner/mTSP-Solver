
#include <iostream>

using namespace std;

class REWARD
{
    // friend class GA;
public:
    REWARD(int timeSteps, double rewardStart, double rewardEnd, double weight)
    {
        _ts = timeSteps;
        _rs = rewardStart;
        _re = rewardEnd;
        _fixedLimit = rewardEnd;
        _w = weight;
    }

    double interactionReward(int step)
    {
        if (step > _rs && step < _fixedLimit)
        {

            double rewards = _w / 2;
            double rewardCoefficient = static_cast<double>(step - _rs) / (_re - _rs);

             if (rewardCoefficient <= 0.5){
            rewards = rewardCoefficient * _w;
             }

            _rs = step;
            _re = _fixedLimit + step;
            return rewards;
        }

        return 0;
    }

    double interactionReward2(int step)
    {
        if (step > _rs && step < _fixedLimit)
        {
            double rewards = _w / 2;
            double rewardCoefficient = static_cast<double>(step - _rs) / (_re - _rs);
             if (rewardCoefficient <= 0.5) {
            rewards = rewardCoefficient * _w;
             }
            return rewards;
        }
        else
        {
            return 0;
        }
    }

private:
    int _ts;
    double _rs, _re, _fixedLimit, _w;
};
