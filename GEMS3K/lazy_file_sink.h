#pragma once
// Rotating file sink that opens its file when the first message reaches it.
//
// spdlog's rotating_file_sink opens (and so creates) its file in the constructor. GEMS3K attaches such a sink to its
// loggers at start-up, so ipmlog.txt and the user's log file would appear even when nothing is ever logged (for
// example with the log level set to off). This sink keeps the same rotation and format but creates the file only
// when there is something to write.

#include <cstddef>
#include <memory>
#include <mutex>
#include <string>
#include <spdlog/sinks/base_sink.h>
#include <spdlog/sinks/rotating_file_sink.h>

namespace gems3k_detail {

template<typename Mutex>
class lazy_rotating_file_sink final : public spdlog::sinks::base_sink<Mutex>
{
public:
    lazy_rotating_file_sink(std::string path, std::size_t max_size, std::size_t max_files)
        : path_(std::move(path)), max_size_(max_size), max_files_(max_files) {}

protected:
    void sink_it_(const spdlog::details::log_msg& msg) override
    {
        open_();
        inner_->log(msg);
    }

    void flush_() override
    {
        if(inner_) inner_->flush();
    }

    // set_pattern() reaches here too: keep the formatter, and pass a copy to the real sink once it exists
    void set_formatter_(std::unique_ptr<spdlog::formatter> sink_formatter) override
    {
        spdlog::sinks::base_sink<Mutex>::set_formatter_(std::move(sink_formatter));
        if(inner_) inner_->set_formatter(this->formatter_->clone());
    }

private:
    void open_()
    {
        if(inner_) return;
        inner_ = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(path_, max_size_, max_files_);
        inner_->set_formatter(this->formatter_->clone());
    }

    std::string path_;
    std::size_t max_size_;
    std::size_t max_files_;
    std::shared_ptr<spdlog::sinks::rotating_file_sink_mt> inner_;
};

using lazy_rotating_file_sink_mt = lazy_rotating_file_sink<std::mutex>;

/// A rotating file sink that creates `path` on the first message written to it.
inline std::shared_ptr<spdlog::sinks::sink> make_lazy_file_sink(const std::string& path,
                                                                std::size_t max_size, std::size_t max_files)
{
    return std::make_shared<lazy_rotating_file_sink_mt>(path, max_size, max_files);
}

} // namespace gems3k_detail
