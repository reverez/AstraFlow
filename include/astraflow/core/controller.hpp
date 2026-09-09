#pragma once
#include "astraflow/analysis/engineering.hpp"
#include <condition_variable>
#include <mutex>
#include <optional>
#include <thread>
namespace astraflow {
enum class RunState { Idle, Ready, Running, Paused, Converged, Finished, Failed };
const char *state_name(RunState);
struct HistoryRow {
    double iteration;
    std::array<double, 4> residual;
};
struct Snapshot {
    std::shared_ptr<const Mesh> mesh;
    std::vector<double> conservative;
    StepStats stats;
    nlohmann::json engineering;
    std::vector<HistoryRow> history;
};
class Controller {
  public:
    explicit Controller(Config);
    ~Controller();
    void run();
    void pause();
    void single_step();
    void reset(Config);
    RunState state() const;
    std::string error() const;
    std::shared_ptr<const Snapshot> snapshot() const;
    bool wait_for(RunState, std::chrono::milliseconds timeout);
    bool wait_iterations(int, std::chrono::milliseconds timeout);

  private:
    void work(std::stop_token);
    void publish(const Simulation &, std::shared_ptr<const Mesh>, const std::vector<HistoryRow> &,
                 RunState);
    mutable std::mutex mutex_;
    std::condition_variable cv_;
    std::optional<Config> pending_;
    bool running_ = false, single_ = false;
    RunState state_ = RunState::Idle;
    std::string error_;
    std::shared_ptr<const Snapshot> snapshot_;
    std::jthread worker_;
};
} // namespace astraflow
