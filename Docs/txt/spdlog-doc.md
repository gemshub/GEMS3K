
##  using [spdlog](https://github.com/gabime/spdlog/tree/master) 

Very fast, header only, C++ logging library. 

source https://github.com/gabime/spdlog/tree/master

documentation https://spdlog.docsforge.com/master/

### Install spdlog

* Ubuntu: apt-get install libspdlog-dev

* Homebrew: brew install spdlog

* Building spdlog library

```sh

        mkdir -p ~/code && \
                cd ~/code && \
                git clone https://github.com/gabime/spdlog -b v1.10.0
        cd spdlog && mkdir -p build && \
                cd build && \
                cmake .. -DCMAKE_CXX_FLAGS=-fPIC && \
                make
        sudo make install

```


### Brief

spdlog levels :  trace = 0, debug = 1, info = 2, warn = 3, err = 4, critical = 5, off = 6



```c++

inline spdlog::level::level_enum from_str(const std::string &name)
{
    static std::unordered_map<std::string, level_enum> name_to_level = // map string->level
        {{level_names[0], level::trace},                               // trace
            {level_names[1], level::debug},                            // debug
            {level_names[2], level::info},                             // info
            {level_names[3], level::warn},                             // warn
            {level_names[4], level::err},                              // err
            {level_names[5], level::critical},                         // critical
            {level_names[6], level::off}};                             // off

    auto lvl_it = name_to_level.find(name);
    return lvl_it != name_to_level.end() ? lvl_it->second : level::off;
}

```

### Used loggers into code

 Thread-safe logger to stdout with colors and/or to file loggers
 
 Currently, such loggers are implemented for gems3k and additional for gems3gui, thermofun, and chemicalfun. It is possible to find loggers by name and change the level or remove logging. 



```c++
std::shared_ptr<spdlog::logger> gems_logger = spdlog::stdout_color_mt("gems3k");

std::shared_ptr<spdlog::logger> TMulti::ipm_logger = spdlog::stdout_color_mt("ipm");

std::shared_ptr<spdlog::logger> TNode::node_logger = spdlog::stdout_color_mt("tnode");

std::shared_ptr<spdlog::logger> TKinMet::kinmet_logger = spdlog::stdout_color_mt("kinmet");

std::shared_ptr<spdlog::logger> TSolMod::solmod_logger = spdlog::stdout_color_mt("solmod");

std::shared_ptr<spdlog::logger> gui_logger = spdlog::stdout_color_mt("gems3gui");

std::shared_ptr<spdlog::logger> thfun_logger = spdlog::stdout_color_mt("thermofun");

std::shared_ptr<spdlog::logger> chfun_logger = spdlog::stdout_color_mt("chemicalfun");
```




 Only file logger ( old ipmlog.txt file)


```c++
std::shared_ptr<spdlog::logger> TNode::ipmlog_file = spdlog::rotating_logger_mt("ipmlog", "ipmlog.txt", 1048576, 3);
```


### API for changing loggers settings in source code


```c++

/// Remove logging to stdout, logging data only to text file logfile_name
    void gems3k_clear_loggers(const std::string &logfile_name);

/// Show/hide logging to stdout (use_cout), to file (hide if logfile_name empty), and set up logging level for all loggers.
    void gems3k_update_loggers(bool use_stdout, const std::string &logfile_name, size_t log_level);

```


Example for gems3kgui


```c++

    gemsSettings().gems3k_update_loggers( true, "gems3k_gui.log", spdlog::level::info);
    gui_logger->set_level(spdlog::level::info);

```

### Using configuration file (*gems3k-config.json*) for changing loggers settings without rebuilding

To read config data add access to ```gemsSettings()``` object.


* **log.level**     - level for all loggers except defined into log.module_level 
* **log.module_level**  - specific level for particular logger
* **log.modules**   - list of existing loggers
* **log.file**  - description of rotating file sink based on size 
* **log.file.modules**  - list of loggers to file
* **log.file.pattern** - pattern for loggers to file (by default *"[%n] [%^%l%$] %v"*)

* **log.stdout.modules** - list of loggers to stdout with colors (by default all)
* **log.stdout.pattern** - pattern for loggers to stdout with colors (by default *"[%n] [%^%l%$] %v"*)



```json

{
    "log": {
        "level": "info",
        "modules": [
            "gems3k",
            "ipm",
            "tnode",
            "kinmet",
            "solmod",
            "thermofun",
            "chemicalfun"
        ],
        "module_level": {
            "chemicalfun": "trace",
            "gems3k": "debug",
            "ipm": "info",
            "thermofun": "trace"
        },
        "file": {
            "count": 10,
            "modules": [
                "ipm",
                "thermofun",
                "chemicalfun"
            ],
            "path": "gems3k-log.txt",
            "size": 10000000
        }
     }
}


```

