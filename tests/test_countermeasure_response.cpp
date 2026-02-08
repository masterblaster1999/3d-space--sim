#include "stellar/sim/MissileDefense.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_countermeasure_response() {
  int fails = 0;

  const math::Vec3d tgtPos{0, 0, 10};
  const math::Vec3d tgtVel{0, 0, 0};

  // --- Discrete seeker updates: align to the next scan tick ---
  {
    sim::Missile m{};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 5};
    m.ttlSimSec = 10.0;
    m.ageSimSec = 0.75;
    m.seeker = sim::MissileSeekerType::Radar;

    m.seekerActivationSimSec = 0.0;
    m.seekerUpdatePeriodSimSec = 0.5;

    sim::CountermeasureResponseParams p{};
    p.decoyWindowSec = 6.0;
    p.releaseLeadSec = 0.05;
    p.panicTtiSec = 0.0; // disable panic so timing behavior is testable
    p.enableRepeatBursts = true;
    p.maxBursts = 3;

    const sim::CountermeasureResponsePlan plan = sim::planCountermeasureResponse(m, tgtPos, tgtVel, p);
    if (!plan.valid) {
      std::cerr << "[test_countermeasure_response] expected discrete plan to be valid.\n";
      ++fails;
    } else {
      if (plan.type != sim::CountermeasureType::Chaff) {
        std::cerr << "[test_countermeasure_response] expected radar seeker -> Chaff recommendation.\n";
        ++fails;
      }

      // age=0.75, period=0.5 => next tick in 0.25s, releaseLead=0.05 => delay ~0.20s.
      if (!approx(plan.firstReleaseDelaySec, 0.20, 1e-6)) {
        std::cerr << "[test_countermeasure_response] unexpected release delay. got=" << plan.firstReleaseDelaySec
                  << " expected ~0.20\n";
        ++fails;
      }

      // Within tti=2, there are 4 scan ticks remaining; maxBursts clamps to 3.
      if (plan.burstCount != 3 || !approx(plan.repeatEverySec, 0.5, 1e-12)) {
        std::cerr << "[test_countermeasure_response] unexpected repeat plan. bursts=" << plan.burstCount
                  << " every=" << plan.repeatEverySec << "\n";
        ++fails;
      }
    }
  }

  // --- Continuous seeker: recommend immediate deployment (no scan alignment) ---
  {
    sim::Missile m{};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 5};
    m.ttlSimSec = 10.0;
    m.ageSimSec = 0.75;
    m.seeker = sim::MissileSeekerType::Heat;

    m.seekerActivationSimSec = 0.0;
    m.seekerUpdatePeriodSimSec = 0.0;

    sim::CountermeasureResponseParams p{};
    p.decoyWindowSec = 6.0;
    p.releaseLeadSec = 0.05;
    p.panicTtiSec = 0.0;
    p.enableRepeatBursts = true;
    p.maxBursts = 3;

    const sim::CountermeasureResponsePlan plan = sim::planCountermeasureResponse(m, tgtPos, tgtVel, p);
    if (!plan.valid) {
      std::cerr << "[test_countermeasure_response] expected continuous plan to be valid.\n";
      ++fails;
    } else {
      if (plan.type != sim::CountermeasureType::Flare) {
        std::cerr << "[test_countermeasure_response] expected heat seeker -> Flare recommendation.\n";
        ++fails;
      }
      if (!approx(plan.firstReleaseDelaySec, 0.0, 1e-12)) {
        std::cerr << "[test_countermeasure_response] expected immediate release for continuous seeker. got="
                  << plan.firstReleaseDelaySec << "\n";
        ++fails;
      }
      if (plan.burstCount != 1 || !approx(plan.repeatEverySec, 0.0, 1e-12)) {
        std::cerr << "[test_countermeasure_response] continuous seeker should not recommend repeats by default. bursts="
                  << plan.burstCount << " every=" << plan.repeatEverySec << "\n";
        ++fails;
      }
    }
  }

  // --- Seeker activation delay: first burst should be scheduled near activation tick ---
  {
    sim::Missile m{};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 5};
    m.ttlSimSec = 10.0;
    m.ageSimSec = 0.2;
    m.seeker = sim::MissileSeekerType::Heat;

    m.seekerActivationSimSec = 1.0;     // seeker comes alive in 0.8s
    m.seekerUpdatePeriodSimSec = 0.5;   // discrete scans

    sim::CountermeasureResponseParams p{};
    p.decoyWindowSec = 6.0;
    p.releaseLeadSec = 0.05;
    p.panicTtiSec = 0.0;
    p.enableRepeatBursts = false;

    const sim::CountermeasureResponsePlan plan = sim::planCountermeasureResponse(m, tgtPos, tgtVel, p);
    if (!plan.valid) {
      std::cerr << "[test_countermeasure_response] expected delayed-seeker plan to be valid.\n";
      ++fails;
    } else {
      // First tick at activation (0.8s), so release at 0.8-0.05 = 0.75.
      if (!approx(plan.firstReleaseDelaySec, 0.75, 1e-6)) {
        std::cerr << "[test_countermeasure_response] activation-aligned delay unexpected. got=" << plan.firstReleaseDelaySec
                  << " expected ~0.75\n";
        ++fails;
      }
    }
  }

  // --- If seeker won't activate before impact, plan is invalid ---
  {
    sim::Missile m{};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 5};
    m.ttlSimSec = 10.0;
    m.ageSimSec = 0.0;
    m.seeker = sim::MissileSeekerType::Heat;

    m.seekerActivationSimSec = 3.0; // activates after predicted impact (~2s)
    m.seekerUpdatePeriodSimSec = 0.5;

    const sim::CountermeasureResponsePlan plan = sim::planCountermeasureResponse(m, tgtPos, tgtVel);
    if (plan.valid) {
      std::cerr << "[test_countermeasure_response] expected invalid plan when seeker activates after impact.\n";
      ++fails;
    }
  }

  return fails;
}
