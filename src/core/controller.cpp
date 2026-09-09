#include "astraflow/core/controller.hpp"
#include <stdexcept>
namespace astraflow {
const char *state_name(RunState s) {
    switch (s) {
    case RunState::Idle:
        return "Preparing";
    case RunState::Ready:
        return "Ready";
    case RunState::Running:
        return "Running";
    case RunState::Paused:
        return "Paused";
    case RunState::Converged:
        return "Converged";
    case RunState::Finished:
        return "Limit reached";
    case RunState::Failed:
        return "Failed";
    }
    return "Unknown";
}
Controller::Controller(Config c) {
    c.validate();
    pending_ = std::move(c);
    worker_ = std::jthread([this](std::stop_token token) { work(token); });
}
Controller::~Controller() {
    worker_.request_stop();
    cv_.notify_all();
    if (worker_.joinable())
        worker_.join();
}
void Controller::run() {
    std::lock_guard lock(mutex_);
    if (state_ != RunState::Ready && state_ != RunState::Paused)
        return;
    running_ = true;
    state_ = RunState::Running;
    cv_.notify_all();
}
void Controller::pause() {
    std::lock_guard lock(mutex_);
    running_ = false;
    cv_.notify_all();
}
void Controller::single_step() {
    std::lock_guard lock(mutex_);
    if (state_ != RunState::Ready && state_ != RunState::Paused)
        return;
    single_ = true;
    cv_.notify_all();
}
void Controller::reset(Config c) {
    c.validate();
    std::lock_guard lock(mutex_);
    pending_ = std::move(c);
    running_ = false;
    single_ = false;
    state_ = RunState::Idle;
    snapshot_.reset();
    error_.clear();
    cv_.notify_all();
}
RunState Controller::state() const {
    std::lock_guard lock(mutex_);
    return state_;
}
std::string Controller::error() const {
    std::lock_guard lock(mutex_);
    return error_;
}
std::shared_ptr<const Snapshot> Controller::snapshot() const {
    std::lock_guard lock(mutex_);
    return snapshot_;
}
bool Controller::wait_for(RunState state, std::chrono::milliseconds timeout) {
    std::unique_lock lock(mutex_);
    return cv_.wait_for(lock, timeout, [&] { return state_ == state; });
}
bool Controller::wait_iterations(int n, std::chrono::milliseconds timeout) {
    std::unique_lock lock(mutex_);
    return cv_.wait_for(lock, timeout,
                        [&] { return snapshot_ && snapshot_->stats.iterations >= n; });
}
void Controller::publish(const Simulation &sim, std::shared_ptr<const Mesh> mesh,
                         const std::vector<HistoryRow> &history, RunState state) {
    auto snapshot = std::make_shared<Snapshot>();
    snapshot->gas = sim.config().settings.gas;
    snapshot->mesh = std::move(mesh);
    snapshot->conservative = sim.state();
    snapshot->stats = sim.stats();
    snapshot->engineering =
        engineering(*snapshot->mesh, snapshot->conservative, sim.config().settings.gas,
                    sim.config().settings.back_pressure);
    snapshot->history = history;
    std::lock_guard lock(mutex_);
    if (pending_)
        return;
    snapshot_ = std::move(snapshot);
    state_ = state;
    cv_.notify_all();
}
void Controller::work(std::stop_token stop) {
    std::unique_ptr<Simulation> sim;
    std::shared_ptr<const Mesh> mesh;
    std::vector<HistoryRow> history;
    while (!stop.stop_requested()) {
        std::optional<Config> request;
        bool advance = false;
        {
            std::unique_lock lock(mutex_);
            cv_.wait(lock, [&] {
                return stop.stop_requested() || pending_ || running_ || single_ ||
                       state_ == RunState::Running;
            });
            if (stop.stop_requested())
                break;
            if (pending_) {
                request = std::move(pending_);
                pending_.reset();
                running_ = false;
                single_ = false;
            } else {
                advance = running_ || single_;
                single_ = false;
            }
        }
        try {
            if (request) {
                sim = std::make_unique<Simulation>(*request);
                mesh = std::make_shared<Mesh>(sim->mesh());
                history.clear();
                publish(*sim, mesh, history, RunState::Ready);
                continue;
            }
            if (!sim)
                continue;
            if (advance) {
                sim->step();
                auto stats = sim->stats();
                history.push_back({double(stats.iterations), stats.residual});
                if (history.size() > 10000)
                    history.erase(history.begin(), history.begin() + 1000);
            }
            bool running = false;
            {
                std::lock_guard lock(mutex_);
                if (sim->finished())
                    running_ = false;
                running = running_;
            }
            RunState state =
                sim->finished()
                    ? (sim->termination() == "converged" ? RunState::Converged : RunState::Finished)
                : running ? RunState::Running
                          : RunState::Paused;
            if (!running || sim->stats().iterations % 10 == 0)
                publish(*sim, mesh, history, state);
        } catch (const std::exception &e) {
            std::lock_guard lock(mutex_);
            error_ = e.what();
            state_ = RunState::Failed;
            running_ = false;
            single_ = false;
            cv_.notify_all();
        }
    }
}
} // namespace astraflow
