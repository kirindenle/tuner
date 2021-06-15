// TODO list
// рефакторинг
// сделать DATA_PATH рантаймовым
// шрифт тоже по-хорошему должен быть рантаймовым !!!!!!!
// адаптивный размер текста! красивый!
// отображение данных в любом месте, а не только снизу сверху или посередине

#define DATA_PATH "C:/Users/22599/projects/_game1/data/"

// defer is from https://gist.github.com/andrewrk/ffb272748448174e6cdb4958dae9f3d8
#define CONCAT_INTERNAL(x,y) x##y
#define CONCAT(x,y) CONCAT_INTERNAL(x,y)

template<typename T>
struct ExitScope {
    T lambda;
    ExitScope(T lambda):lambda(lambda){}
    ~ExitScope(){lambda();}
    ExitScope(const ExitScope&);
  private:
    ExitScope& operator =(const ExitScope&);
};

class ExitScopeHelp {
  public:
    template<typename T>
        ExitScope<T> operator+(T t){ return t;}
};

#define defer const auto& CONCAT(defer__, __LINE__) = ExitScopeHelp() + [&]()

// ************************************ START ************************************




#include <assert.h>
#include "array.h"

#include "SDL.h"
#include "SDL_ttf.h"
#include <stdio.h>
#include <thread>
#include <atomic>
#include <string>
#include <algorithm>
#include <cmath>

using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

using u64atom = std::atomic_uint64_t;

static constexpr int DEFAULT_SAMPLES = 4096;
static constexpr int RECORD_BUFFER_SIZE = DEFAULT_SAMPLES * 128;

void set_default_spec(SDL_AudioSpec *spec);

struct Devices {
    char const* play_device   = NULL;
    char const* record_device = NULL;
};
bool get_play_and_record_devices(Devices* devices_to_fill, int n_args, char const* const* args);

struct Record {
    u64atom current_writen_pos = 0;
    Array<u8> record_buffer;
};

struct Audio_Wav {
    Array<u8>     wav;
    SDL_AudioSpec spec;
    char const*   name = NULL;
};

struct Player_User_Data {
    Record* record;
    u64     position = 0;
};

struct Color : public SDL_Color {
    Color (u32 code) {
        a = code & 0xFF; code >>= 8;
        b = code & 0xFF; code >>= 8;
        g = code & 0xFF; code >>= 8;
        r = code & 0xFF;
    }
    Color (u8 _r, u8 _g, u8 _b, u8 _a) {
        r = _r;
        g = _g;
        b = _b;
        a = _a;
    }
};

struct Context;
struct Note {
    char const* letter;
    int         octave;
};

void play_from_mic(Devices* devices);
void visualize_from_mic(Devices* devices);
void visualize_wav(Audio_Wav const* audio);

bool play_wav_from_drive(char const* wav_file, char const* play_device);
void play_audio(Player_User_Data * player_user_data, u8 * stream, int chunk_len);
void record_audio(Record * record, u8 * stream, int chunk_len);
bool open_wav(Audio_Wav * audio_wav, char const* file);
bool open_my_wav(Array<float> * audio, char const* file);
void draw_vertical_line(Context *context, Color c, double pos);
bool testing_weird_zcr(Context * context, Array<float> const& audio);
float Get16bitAudioSample(u8 const* bytebuffer, SDL_AudioFormat format); // TODO переписать это для всех форматов, вроде есть CRV какой-то, а может блин просто надо было изначально float32 просить
bool init_sdl(Context * context);
void close_sdl(Context * context);
double get_spline_value(Array<float> const& audio, double x);
double hz_to_units(double hz);
Note units_to_note(double units);
// dddddddddddddddd

double hz_to_units(double hz) {
    // using log2
    static constexpr double log_of_2   = 1.0;
    static constexpr double log_of_440 = 8.781359713524659938;
    return 12 * (std::log2(hz) - log_of_440) / log_of_2;
}

Note units_to_note(double units) {
    static char const* notes[12] = {"A", "A#", "H", "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#"};
    int closest_note = (int) std::round(units);
    int mod_12;
    if (closest_note < 0) {
        mod_12 = (closest_note + 12 * std::abs(closest_note / 12) + 12) % 12;
    } else {
        mod_12 = closest_note % 12;
    }
    Note ret;
    ret.letter = notes[mod_12];
    ret.octave = (closest_note - mod_12) / 12 + 4; // 440 hz is A4
    return ret;
}

struct Context {
    SDL_Window   *window;
    int win_w, win_h;
    SDL_Renderer *renderer;
    TTF_Font     *font;
    SDL_AudioSpec spec;
};


int main(int const n_args, char ** const args) {
    if(SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO) == -1) {SDL_Log("Could not initialize SDL: %s.\n", SDL_GetError()); return -1;}
    defer {SDL_Quit();};

    Devices devices;
    bool got_devices = get_play_and_record_devices(&devices, n_args, args);
    if (!got_devices) return -1;

    // play_wav_from_drive("440Hz.wav", devices.play_device);
    // play_from_mic(&devices);


    visualize_from_mic(&devices);

    // Audio_Wav audio_wav;
    // if (!open_wav(&audio_wav, "440Hz-20dB.wav")) return -1;
    // defer {SDL_FreeWAV(audio_wav.wav.data);};
    // visualize_wav(&audio_wav);

//    static float storage_my_wav[DEFAULT_SAMPLES];
//    Array<float> my_wav_audio(storage_my_wav, 0, DEFAULT_SAMPLES);
//    open_my_wav(&my_wav_audio, "record25");
//
//    Context storage_context;
//    Context *context = &storage_context;
//    init_sdl(context);
//    defer {close_sdl(context);};
//
//    testing_weird_zcr(context, my_wav_audio);

    return 0;
}

void set_default_spec(SDL_AudioSpec *spec) {
    spec->freq     = 44100;
    spec->channels = 1;
    spec->samples  = DEFAULT_SAMPLES;
    spec->callback = NULL;
    spec->userdata = NULL;
    spec->format   = AUDIO_S16;
}

void clamp(int * var, int min, int max) {
    assert(var);
    if (*var > max) *var = max;
    if (*var < min) *var = min;
}

bool open_wav(Audio_Wav * audio_wav, char const* file) {
    std::string relative_path_to_wav_file = DATA_PATH;
    relative_path_to_wav_file += file;
    {
        u32 size;
        SDL_AudioSpec * tmp = SDL_LoadWAV(relative_path_to_wav_file.c_str(), &audio_wav->spec, &audio_wav->wav.data, &size);
        if(tmp != &audio_wav->spec) {SDL_Log("Could not open audio file: %s\n", relative_path_to_wav_file.c_str()); return false;}
        assert(size <= INT_MAX);
        audio_wav->wav.size     = size;
        audio_wav->wav.capacity = size;
    }
    audio_wav->name = file;

    SDL_AudioSpec const* spec = &audio_wav->spec;
    SDL_Log("Opened file %s in format freq=%d, format=%x, channels=%d, samples=%d", file, spec->freq, spec->format, spec->channels, spec->samples);
    return true;
}

bool open_my_wav(Array<float> * audio, char const* filename) {
    assert(audio);
    assert(filename);
    std::string full_path = DATA_PATH;
    full_path += filename;
    static u8 storage_raw[DEFAULT_SAMPLES * 2];
    int data_size = -1;
    {
        SDL_RWops *file = SDL_RWFromFile(full_path.c_str(), "r");
        if (!file) {SDL_Log("Error opening file to read %s: %s", full_path.c_str(), SDL_GetError()); return false;}
        defer {SDL_RWclose(file);};

        size_t read = SDL_RWread(file, &data_size, sizeof(int), 1);
        assert(data_size <= DEFAULT_SAMPLES * 2);
        assert(data_size >= 0);
        if (read != 0) read = SDL_RWread(file, &storage_raw, data_size, 1);
        if (read == 0) {SDL_Log("Error reading file %s: %s", full_path.c_str(), SDL_GetError()); return false;}
    }
    assert(data_size != -1);

    // чот некраиво
    SDL_AudioSpec default_spec;
    set_default_spec(&default_spec);

    assert(data_size % 2 == 0);
    assert(audio->capacity >= data_size / 2);
    audio->size = data_size / 2;
    for (int i = 0; i < data_size / 2; ++i) {
        (*audio)[i] = Get16bitAudioSample(&storage_raw[2 * i], default_spec.format); // TODO перенести spec в запись
    }
    return true;
}

void record_audio(Record * record, u8 * stream, int chunk_len) {
    // static std::chrono::high_resolution_clock clock;
    // static auto last_sample_moment = clock.now();
    // auto current_moment = clock.now();
    // SDL_Log("%f", std::chrono::duration<double>(current_moment - last_sample_moment).count());
    // last_sample_moment = current_moment;

    u64 start_pos = record->current_writen_pos.load();
    // SDL_Log("record audio %d %d", (int)start_pos, (int)chunk_len);
    if (start_pos + chunk_len > RECORD_BUFFER_SIZE) {
        start_pos = 0;
    }
    memcpy(&record->record_buffer[start_pos], stream, chunk_len);
    record->current_writen_pos.store(start_pos + chunk_len);
}

bool get_play_and_record_devices(Devices* devices_to_fill, int const n_args, char const* const* const args) {
    assert(devices_to_fill);
    if (n_args != 1 && n_args != 3) {SDL_Log("need 0 or 2 arguments: idx_playback idx_record"); return false;}

    int play_device_idx   = -1;
    int record_device_idx = -1;
    int n_play_devices   = SDL_GetNumAudioDevices(SDL_FALSE);
    int n_record_devices = SDL_GetNumAudioDevices(SDL_TRUE);
    if (n_args == 1) {
        play_device_idx   = 0;
        record_device_idx = 0;
        assert(play_device_idx < n_play_devices);
        assert(record_device_idx < n_record_devices);
        // SDL_Log("presented audio playback devices:");
        // for (int i = 0; i < n_play_devices; ++i) {
        //     SDL_Log("%d: %s", i, SDL_GetAudioDeviceName(i, SDL_FALSE));
        // }
        // SDL_Log("presented audio recording devices:");
        // for (int i = 0; i < n_record_devices; ++i) {
        //     SDL_Log("%d: %s", i, SDL_GetAudioDeviceName(i, SDL_TRUE));
        // }
        // return false;
    }
    if (n_args == 3) {
        play_device_idx = SDL_atoi(args[1]);
        if (play_device_idx == 0 && args[1][0] != '0'
          ||play_device_idx >= n_play_devices) {SDL_Log("you should choose playback device number from 0 to %d, but you entered %s", n_play_devices - 1, args[1]); return false;}
        record_device_idx = SDL_atoi(args[2]);
        if (record_device_idx == 0 && args[2][0] != '0'
          ||record_device_idx >= n_record_devices) {SDL_Log("you should choose recording device number from 0 to %d, but you entered %s", n_record_devices - 1, args[2]); return false;}
    }
    if (play_device_idx == -1)   {SDL_Log("Pass the correct number as second argument but not '%s'", args[1]); return false;}
    if (record_device_idx == -1) {SDL_Log("Pass the correct number as third  argument but not '%s'", args[2]); return false;}

    devices_to_fill->play_device   = SDL_GetAudioDeviceName(play_device_idx, SDL_FALSE);
    devices_to_fill->record_device = SDL_GetAudioDeviceName(record_device_idx, SDL_TRUE);

    SDL_Log("You chose playback device %s",  devices_to_fill->play_device);
    SDL_Log("You chose recording device %s", devices_to_fill->record_device);
    return true;
}

// *************************** PLAY WAV STUFF ***************************
struct Player_State {
    u32  rest_length;
    u8 * play_position;
};

void play_audio_data(void * userdata, u8 * stream, int chunk_len) {
    assert(userdata);
    Player_State * audio_data = (Player_State *) userdata;
    if (audio_data->rest_length == 0) return;

    clamp(&chunk_len, 0, audio_data->rest_length);
    SDL_memcpy (stream, audio_data->play_position, chunk_len); // simply copy from one buffer into the other
    audio_data->play_position += chunk_len;
    audio_data->rest_length   -= chunk_len;
}


// ****************************************** visualization ******************************************

// from https://github.com/arcogabbo/musicvisualizer/blob/master/data.c
float Get16bitAudioSample(u8 const* bytebuffer, SDL_AudioFormat format)
{
    u16 val= 0;

    if(SDL_AUDIO_ISLITTLEENDIAN(format)) val = (uint16_t)bytebuffer[0] | ((uint16_t)bytebuffer[1] << 8);
    else                                 val = ((uint16_t)bytebuffer[0] << 8) | (uint16_t)bytebuffer[1];

    if(SDL_AUDIO_ISSIGNED(format)) return ((int16_t)val)/32768.0f;
    else                           return val/65535.0f;
}

enum Graph_Pos {
    TOP,
    BOTTOM,
    MIDDLE
};

struct Sorted_Peaks;
struct Context;

void data_to_points(SDL_Point* points /* size is window_w */, const float* data, int size, int window_w, int window_h, Graph_Pos pos);
void calculate_peaks(Sorted_Peaks * peaks, float * data, float * data_img, int n_data);
// void bubble_sort(int * ids, float * vals, int size);
void draw_text(Context *context, char const* text, int x, int y);

static unsigned print_slower = 0;


void find_frequency_by_fft(Context*, const Array<float>&);
void find_frequency_by_zcr(Context*, const Array<float>&, bool print_intervals = true);
bool wait_until_press(SDL_Keycode * pressed_key = NULL);
bool selection_mode(Context * context, Array<float> * audio, Audio_Wav * wav_raw = NULL);
void draw_audio(Context * context, Array<float> const& audio);

void visualize_from_mic(Devices* devices) {
    Context storage_context;
    Context * context = &storage_context;
    if (!init_sdl(context)) return;
    defer {close_sdl(context);};

    set_default_spec(&context->spec);
    context->spec.callback = (SDL_AudioCallback)record_audio;

    static u8 record_buffer[RECORD_BUFFER_SIZE];
    Record record;
    record.record_buffer.data = record_buffer;
    record.record_buffer.size = RECORD_BUFFER_SIZE;
    record.record_buffer.capacity = RECORD_BUFFER_SIZE;

    context->spec.userdata = &record; // TODO typed callbacks

    SDL_AudioDeviceID record_id;
    record_id = SDL_OpenAudioDevice(devices->record_device, SDL_TRUE, &context->spec, NULL, 0);
    if (record_id == 0) {SDL_Log("Failed to open audio record device: %s", SDL_GetError()); return; } defer {SDL_CloseAudioDevice(record_id);};

    SDL_PauseAudioDevice(record_id, 0);

    double dt = 0;
    std::chrono::high_resolution_clock clock;
    auto prev_time = clock.now();

    const double sample_time = (double) context->spec.samples / context->spec.freq;
    auto last_sample_moment = clock.now();
    // TODO синхронизация с мирофоном
    for(;;) {
        // SDL_Log("sample_tome %f", sample_time);
        auto current_time = clock.now();
        dt = std::chrono::duration<double>(current_time - prev_time).count();
        prev_time = current_time;


        double time_since_last_sample = std::chrono::duration<double>(current_time - last_sample_moment).count();
        if (time_since_last_sample >= sample_time) {
            //SDL_Log("time_since_last_sample = %f", time_since_last_sample);
            last_sample_moment = current_time;
        }

        static u8 storage_samples_raw [DEFAULT_SAMPLES * 2];
        Audio_Wav wav_raw;
        wav_raw.spec = context->spec; // TODO убрать spec из контекста
        // ебаные конструкторы, просто так не вызовешь на заранее объявленной локальной переменной
        wav_raw.wav.data = storage_samples_raw;
        wav_raw.wav.capacity = wav_raw.wav.size = DEFAULT_SAMPLES * 2;
        {
            // TODO компилятор тут не ругается, что использована неинициализированная переменная, кек
            u64 record_end_pos = record.current_writen_pos.load();
            while (record_end_pos < DEFAULT_SAMPLES * 2) {
                record_end_pos = record.current_writen_pos.load();
            }
            u64 record_read_pos = record_end_pos - DEFAULT_SAMPLES * 2;
            SDL_memcpy(storage_samples_raw, record_buffer + record_read_pos, DEFAULT_SAMPLES * 2);
        }

        static float storage_samples [DEFAULT_SAMPLES];
        Array<float> audio(storage_samples, DEFAULT_SAMPLES);

        for (int i = 0; i < DEFAULT_SAMPLES; ++i) {
            audio[i] = Get16bitAudioSample(&storage_samples_raw[2 * i], context->spec.format); // TODO перенести spec в запись
        }

        // SDL_Log("dt = %f", dt);
        ++print_slower;

        SDL_SetRenderDrawColor(context->renderer, 0, 0, 0, 255);
        SDL_RenderClear(context->renderer);

        draw_audio(context, audio);

        // find_frequency_by_fft(context, audio);
        find_frequency_by_zcr(context, audio);

        SDL_RenderPresent(context->renderer);
        SDL_Event event;
        bool quit = false;
        while (SDL_PollEvent(&event)) {
            if      (event.type == SDL_QUIT) quit = true;
            else if (event.type == SDL_KEYUP) {
                // TODO SDL_KeyCode и SDL_Keycode можно ли сравнивать?
                if (event.key.keysym.sym == SDLK_RETURN) {
                    // TODO просто сохранять выбранное в отдельный буфер без использования файлов
                    //      и переходить в состояние, в котором вместо чтения с микрофона я буду читать из этого буфера
                    quit = selection_mode(context, &audio, &wav_raw);
                    // TODO избавиться от этих констант, наверное, или назвать как-то нормально
                    if (!quit) {
                        if (audio.data != storage_samples || audio.size != DEFAULT_SAMPLES) {
                            quit = testing_weird_zcr(context, audio);
                        }
                    }
                } else {
                    quit = wait_until_press();
                    break;
                }
            }
        }
        if (quit) return;
    }
}

void testing_clear_screen(Context * context, Array<float> const& audio) {
    SDL_SetRenderDrawColor(context->renderer, 0, 0, 0, 255); SDL_RenderClear(context->renderer);
    find_frequency_by_zcr(context, audio, false);
    draw_audio(context, audio);
    SDL_RenderPresent(context->renderer);
}

float global_y_scale = 1.0f;

void testing_find_frequency_by_zcr(Context * context, Array<float> const& audio, bool print_intervals) {
    static int    intervals_appearences[DEFAULT_SAMPLES];
    static double intervals_sum_diff   [DEFAULT_SAMPLES];
    char text_buffer[100];
    for (int i = 0; i < DEFAULT_SAMPLES; ++i) intervals_appearences[i] = 0;
    for (int i = 0; i < DEFAULT_SAMPLES; ++i) intervals_sum_diff   [i] = 0;
    // memset(intervals_appearences, 0, sizeof(int) * DEFAULT_SAMPLES);

    int    prev_zero       = 0;
    double prev_zero_exact = 0.0;
    // bool first = true;

    For(audio) {
        if (it_index + 1 == audio.size) break;

        auto it_next = audio[it_index + 1];
        if (!(it < 0.0f && it_next >= 0.0f)) continue;

        // SDL_Log("%f %f", audio[i], audio[i + 1]);
        double current_zero_exact = it;
        if (it_next - it >= -it) {
            current_zero_exact = it_index + (-it / (it_next - it));
        } else continue;

        if (prev_zero != 0) {
            const int interval = it_index - prev_zero;
            const double interval_exact = current_zero_exact - prev_zero_exact;

            // отфильтровать интервал, проверить периодичность
            bool bad_zero_cross = false;
            {
                using std::fabs;
                float amplitude = 0.0f;
                int amplitude_idx = 0;
                FOR(range_idx, 0, interval) {
                    if (it_index + range_idx >= audio.size) break;
                    float val = fabs(audio[it_index + range_idx]);
                    if (val > amplitude) {
                        amplitude = val;
                        amplitude_idx = it_index + range_idx;
                    }
                }

                const float tolerance = 0.1f * amplitude;

                double step = interval / interval_exact;
                FOR(range_idx, 0, interval) {
                    if (it_index + range_idx >= audio.size) break;
                    double x_left  = prev_zero_exact + range_idx * step;
                    double x_right = x_left + (current_zero_exact - prev_zero_exact);

                    double y_left  = get_spline_value(audio, x_left);
                    double y_right = get_spline_value(audio, x_right);

                    for (;;) {
                        testing_clear_screen(context, audio);
                        draw_vertical_line(context, 0xFFFFFFFF, float(amplitude_idx) / audio.size);
                        draw_vertical_line(context, 0x00FF00FF, x_left / audio.size);
                        draw_vertical_line(context, 0x0000FFFF, x_right / audio.size);

                        snprintf(text_buffer, 100, "f(%.1e) = %.1e", x_left, y_left);          draw_text(context, text_buffer, 0, 0);
                        snprintf(text_buffer, 100, "f(%.1e) = %.1e | %.1e | %.1e", x_right, y_right, fabs(y_left - y_right), fabs(y_left - y_right) - tolerance);        draw_text(context, text_buffer, 0, 36);
                        snprintf(text_buffer, 100, "amp = %.1e, int = %d, step = %f", amplitude, interval, step); draw_text(context, text_buffer, 0, 72);
                        snprintf(text_buffer, 100, "prev from %f to %f", prev_zero_exact, prev_zero_exact + step * (interval-1)); draw_text(context, text_buffer, 0, 108);
                        snprintf(text_buffer, 100, "curr from %f to %f", current_zero_exact, current_zero_exact + step * (interval-1)); draw_text(context, text_buffer, 0, 144);
                        snprintf(text_buffer, 100, "prev_zero = %d, prev_zero_exact = %f", prev_zero, prev_zero_exact); draw_text(context, text_buffer, 0, 180);
                        snprintf(text_buffer, 100, "current_zero = %d, current_zero_exact = %f", it_index, current_zero_exact); draw_text(context, text_buffer, 0, 216);

                        SDL_RenderPresent(context->renderer);
                        SDL_Keycode key;
                        if (wait_until_press(&key)) exit(0);
                        if (key == SDLK_UP) {
                            global_y_scale *= 1.1f;
                            continue;
                        }
                        if (key == SDLK_DOWN) {
                            global_y_scale /= 1.1f;
                            continue;
                        }
                        break;
                    }

                    if (fabs(y_left - y_right) > tolerance) {
                        bad_zero_cross = true;
                        SDL_Log("bad zero: %d", it_index);
                        draw_vertical_line(context, 0xFF0000FF, current_zero_exact / audio.size);
                        // if (first) {
                        //     draw_vertical_line(context, 0xFFFFFFFF, x_left  / audio.size);
                        //     draw_vertical_line(context, 0xFFFFFFFF, x_right / audio.size);
                        //     first = false;
                        // }
                        goto bad_zero;
                    }
                }
            }
            if (bad_zero_cross) continue;


            snprintf(text_buffer, 100, "adding intexrval %d", interval);   draw_text(context, text_buffer, 0, context->win_h/2 - 18);
            SDL_RenderPresent(context->renderer);
            if (wait_until_press()) exit(0);

            ++intervals_appearences[interval];
            intervals_sum_diff[interval] += (current_zero_exact - prev_zero_exact);
        }

        draw_vertical_line(context, 0x804040FF, (float) it_index / audio.size);
        prev_zero = it_index;
        prev_zero_exact = current_zero_exact;

        bad_zero: ;
    }

    int most_often_interval = 0;
    for (int i = 0; i < DEFAULT_SAMPLES; ++i) {
        if (intervals_appearences[i] > intervals_appearences[most_often_interval]) most_often_interval = i;
    }
    double most_often_interval_exact = intervals_sum_diff[most_often_interval] / intervals_appearences[most_often_interval]; // take avg

    double storage_intervals[DEFAULT_SAMPLES];
    Array <double> intervals(storage_intervals, 0, DEFAULT_SAMPLES);

    int storage_intervals_ids[DEFAULT_SAMPLES];
    Array <int> intervals_ids(storage_intervals_ids, 0, DEFAULT_SAMPLES);

    for (int i = 0; i < DEFAULT_SAMPLES; ++i) {
        if (intervals_appearences[i] > 0) {
            intervals.add(intervals_sum_diff[i] / intervals_appearences[i]);
            intervals_ids.add(i);
        }
    }

    static char buffer[100];
    int last_y = 0;
    if (print_intervals) {
        For (intervals) {
            SDL_Point line[2];
            line[0].x = (int) (it / DEFAULT_SAMPLES * context->win_w);
            line[0].y = 0;
            line[1].x = line[0].x;
            line[1].y = context->win_h;

            // TODO(den): рисовать тут горизонтальные линии красивые
            // SDL_SetRenderDrawColor(context->renderer, 0, 255, 0, 255); SDL_RenderDrawLines(context->renderer, line, 2);
            double Hz = context->spec.freq / it;
            SDL_snprintf(buffer, 100, "%d) %1.2f     %d", intervals_ids[it_index], Hz, intervals_appearences[intervals_ids[it_index]]);
            draw_text(context, buffer, line[0].x, last_y);

            last_y += 30;
        }
    }

    // SDL_Log("%d", most_often_interval);
    // double Hz = 1.0 / ((double) most_often_interval / context->spec.freq);
    double Hz = context->spec.freq / most_often_interval_exact;
    static char buffer_res[100];
    if (buffer_res[0] == 0 || Hz < 1200) {
        SDL_snprintf(buffer_res, 100, "%1.2f", Hz);
    }
    draw_text(context, buffer_res, context->win_w - 200, context->win_h - 60);
}

bool testing_weird_zcr(Context * context, Array<float> const& audio) {
    for(;;) {
        testing_clear_screen(context, audio);
        testing_find_frequency_by_zcr(context, audio, false);
        SDL_RenderPresent(context->renderer);

        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if      (event.type == SDL_QUIT) return true;
            else if (event.type == SDL_KEYUP) {
                if (event.key.keysym.sym == SDLK_ESCAPE) {
                    return false;
                }
            }
        }
    }
}

bool init_sdl(Context * context) {
    bool success = false;
    context->window  = SDL_CreateWindow("Audio visualizer", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 1024, 768, SDL_WINDOW_OPENGL);
    if (!context->window) {SDL_Log("Could not create window: %s\n", SDL_GetError()); return false;}
    SDL_GetWindowSize(context->window, &context->win_w, &context->win_h);
    defer {if (!success) SDL_DestroyWindow(context->window);};

    context->renderer = SDL_CreateRenderer(context->window, -1, SDL_RENDERER_ACCELERATED);
    if (!context->renderer) {SDL_Log("Could not create rendered: %s\n", SDL_GetError()); return false;}
    SDL_SetRenderDrawBlendMode(context->renderer, SDL_BLENDMODE_BLEND);
    defer {if (!success) SDL_DestroyRenderer(context->renderer);};

    if (TTF_Init() == -1) {SDL_Log("Could not init TTF: %s\n", SDL_GetError()); return false;}
    defer {if (!success) TTF_Quit();};

    char const *font_path = DATA_PATH "Poppins-Regular.ttf";
    context->font = TTF_OpenFont(font_path, 36);
    if (!context->font) {SDL_Log("Can't open font %s: %s", font_path, SDL_GetError()); return false;}
    defer {if (!success) TTF_CloseFont(context->font);};

    success = true;
    return success;
}

void close_sdl(Context * context) {
    defer {SDL_DestroyWindow(context->window);};
    defer {SDL_DestroyRenderer(context->renderer);};
    defer {TTF_Quit();};
    defer {TTF_CloseFont(context->font);};
}

bool wait_until_press(SDL_Keycode * pressed_key) {
    SDL_Event event;
    for (;;) {
        if (SDL_PollEvent(&event)) {
            // TODO сделать версию с keydown
            if (event.type == SDL_KEYUP) {
                if (pressed_key) *pressed_key = event.key.keysym.sym;
                return false;
            }
            if (event.type == SDL_QUIT)  return true;
        }
    }
}

// size bytes
// [int][data]
std::string save_my_wav(Audio_Wav const& wav_raw, int start_idx, int end_idx) {
    #define SAVE_RECORD_NAME "record"
    SDL_RWops *file = NULL;
    std::string filename = "";
    for (int i = 0; i < 100; ++i) {
        if (i < 10) filename = DATA_PATH SAVE_RECORD_NAME "0" + std::to_string(i);
        else        filename = DATA_PATH SAVE_RECORD_NAME     + std::to_string(i);

        file = SDL_RWFromFile(filename.c_str(), "r");
        if (file) {
            SDL_RWclose(file);
            continue;
        }

        file = SDL_RWFromFile(filename.c_str(), "wb");
        defer {SDL_RWclose(file);};
        if (!file) {SDL_Log("Error opening file to write %s: %s", filename.c_str(), SDL_GetError()); return "";}

        size_t writen = SDL_RWwrite(file, &wav_raw.wav.size, sizeof(int), 1);
        if (writen != 0) writen = SDL_RWwrite(file, wav_raw.wav.data, wav_raw.wav.size, 1);
        if (writen == 0) {SDL_Log("Error writing to file %s: %s", filename.c_str(), SDL_GetError()); return "";}

        return filename;
    }

    SDL_Log("too much records in folder %s", DATA_PATH);
    return "";
    #undef SAVE_RECORD_NAME
}

bool selection_mode(Context * context, Array<float> * audio, Audio_Wav * wav_raw) {
    int selected_x[2];
    selected_x[0] = 0;
    selected_x[1] = 0;

    bool drag = false;
    int  last_drag_idx = -1;
    bool something_changed = true;
    std::string saved_filename = "";
    for (;;) {
        SDL_Event event;
        if (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT)  return true;
            if (event.type == SDL_KEYUP) {
                if (event.key.keysym.sym == SDLK_RETURN) return false;
                if (event.key.keysym.sym == 's') {
                    int start_audio_idx = int((float) selected_x[0] / context->win_w * audio->size);
                    int end_audio_idx   = int((float) selected_x[1] / context->win_w * audio->size);

                    if (wav_raw) {
                        int raw_bytesize = SDL_AUDIO_BITSIZE(wav_raw->spec.format) / 8;
                        int start_idx = start_audio_idx * raw_bytesize;
                        int end_idx   = end_audio_idx   * raw_bytesize;

                        saved_filename = save_my_wav(*wav_raw, start_idx, end_idx);

                        draw_text(context, "saved record to", 0, 0);
                        draw_text(context, saved_filename.c_str() , 0, 30);
                        SDL_RenderPresent(context->renderer);

                        SDL_Keycode pressed_key;
                        bool quit = wait_until_press(&pressed_key);
                        if (quit) return true;
                        if (pressed_key != SDLK_ESCAPE) {
                            audio->data += start_audio_idx;
                            audio->size = end_audio_idx - start_audio_idx;
                            audio->capacity = end_audio_idx - start_audio_idx;

                            return false;
                        }

                        // TODO может не надо тут?
                        something_changed = true;
                    }

                }
            }

            if (event.type & SDL_MOUSEMOTION) {
                if (event.type == SDL_MOUSEBUTTONDOWN) {
                    drag = true;
                }

                if (event.type == SDL_MOUSEBUTTONUP) {
                    drag = false;
                    last_drag_idx = -1;
                }

                if (!drag) continue;
                something_changed = true;

                int drag_idx;
                int x = event.motion.x;
                clamp(&x, 0, context->win_w - 1);
                if (x > selected_x[1]) drag_idx = 1;
                if (x < selected_x[0]) drag_idx = 0;
                if (x - selected_x[0] < selected_x[1] - x) drag_idx = 0;
                else                                   drag_idx = 1;

                if (last_drag_idx == -1) last_drag_idx = drag_idx;
                if (last_drag_idx != drag_idx) selected_x[last_drag_idx] = selected_x[drag_idx];

                selected_x[drag_idx] = x;
                last_drag_idx = drag_idx;
            }
        }

        while (something_changed) {
            something_changed = false;
            SDL_SetRenderDrawColor(context->renderer, 0, 0, 0, 255);
            SDL_RenderClear(context->renderer);

            find_frequency_by_zcr(context, *audio, false);
            draw_audio(context, *audio);

            SDL_Rect selected_rect;
            selected_rect.x = selected_x[0];
            selected_rect.y = 0;
            selected_rect.w = selected_x[1] - selected_x[0];
            selected_rect.h = context->win_h;
            SDL_SetRenderDrawColor(context->renderer, 255, 0, 255, 64); SDL_RenderFillRect(context->renderer, &selected_rect);
            draw_vertical_line(context, 0xFF00FFFF, (float) selected_x[0] / context->win_w);
            draw_vertical_line(context, 0xFF00FFFF, (float) selected_x[1] / context->win_w);

            SDL_RenderPresent(context->renderer);
        }
    }
}

void draw_audio(Context * context, Array<float> const& audio) {
    static SDL_Point points_orig [DEFAULT_SAMPLES]; // size = width просто я хочу сделать resize поэтому не выношу width в константы
    float const* samples_orig = audio.data;
    data_to_points(points_orig, samples_orig, audio.size, context->win_w, context->win_h, MIDDLE);
    SDL_SetRenderDrawColor(context->renderer, 0, 0, 255, 255); SDL_RenderDrawLines(context->renderer, points_orig, context->win_w);
    SDL_SetRenderDrawColor(context->renderer, 0, 255, 0, 255); SDL_RenderDrawLine(context->renderer, 0, context->win_h/2, context->win_w, context->win_h/2);
}

void draw_text(Context *context, char const* text, int x, int y) {//, char const* text, int x, int y) {
    SDL_Rect message_rect;
    message_rect.x = x + 2;
    message_rect.y = y + 2;
    TTF_SizeText(context->font, text, &message_rect.w, &message_rect.h);

    {
        SDL_Surface* surfaceMessage = TTF_RenderText_Solid(context->font, text, {128, 128, 128});
        if (!surfaceMessage) {SDL_Log("error: %s", SDL_GetError()); return;} defer {SDL_FreeSurface(surfaceMessage);};

        SDL_Texture* texture = SDL_CreateTextureFromSurface(context->renderer, surfaceMessage);
        if (!texture) {SDL_Log("SDL_CreateTextureFromSurface: %s", SDL_GetError()); return;} defer {SDL_DestroyTexture(texture);};

        SDL_RenderCopy(context->renderer, texture, NULL, &message_rect);
    }

    message_rect.x = x;
    message_rect.y = y;

    {
        SDL_Surface* surfaceMessage = TTF_RenderText_Solid(context->font, text, {255, 255, 255});
        if (!surfaceMessage) {SDL_Log("error: %s", SDL_GetError()); return;} defer {SDL_FreeSurface(surfaceMessage);};

        SDL_Texture* texture = SDL_CreateTextureFromSurface(context->renderer, surfaceMessage);
        if (!texture) {SDL_Log("SDL_CreateTextureFromSurface: %s", SDL_GetError()); return;} defer {SDL_DestroyTexture(texture);};

        SDL_RenderCopy(context->renderer, texture, NULL, &message_rect);
    }
}

void data_to_points(SDL_Point* points /* size is window_w */, const float* data, int size, int window_w, int window_h, Graph_Pos pos) {
    bool abs_and_scale = true;
    if (pos == MIDDLE) abs_and_scale = false;
    float scale = 200.0f;
    for (int x = 0; x < window_w; ++x) {
        points[x].x = x;

        int idx = int((float) x / window_w * size);
        float h = data[idx] * global_y_scale;
        if (abs_and_scale) h /= scale;
                if (pos == TOP) points[x].y = int(h * window_h);
        else if (pos == BOTTOM) points[x].y = int(window_h - h * window_h);
        else if (pos == MIDDLE) points[x].y = int((window_h/2) - h * (window_h/2));
        else points[x].y = window_h - 10;
    }
}


double get_spline_value(Array<float> const& audio, double x) {
    assert(x >= 0.0f);
    int    i = (int) x; // floor
    double e = x - i;
    return audio[i] + (audio[i+1] - audio[i]) * e;
}

void draw_vertical_line(Context *context, Color c, double pos) {
    assert(pos >= 0);
    assert(pos <= 1);
    SDL_SetRenderDrawColor(context->renderer, c.r, c.g, c.b, c.a);
    int screen_pos = (int) (pos * context->win_w);
    SDL_RenderDrawLine(context->renderer, screen_pos, 0, screen_pos, context->win_h);
}

void find_frequency_by_zcr(Context * context, Array<float> const& audio, bool print_intervals) {
    static int    intervals_appearences[DEFAULT_SAMPLES];
    static double intervals_sum_diff   [DEFAULT_SAMPLES];
    for (int i = 0; i < DEFAULT_SAMPLES; ++i) intervals_appearences[i] = 0;
    for (int i = 0; i < DEFAULT_SAMPLES; ++i) intervals_sum_diff   [i] = 0;
    // memset(intervals_appearences, 0, sizeof(int) * DEFAULT_SAMPLES);

    int    prev_zero       = 0;
    double prev_zero_exact = 0.0;
    // bool first = true;

    For(audio) {
        if (it_index + 1 == audio.size) break;

        auto it_next = audio[it_index + 1];
        if (!(it < 0.0f && it_next >= 0.0f)) continue;

        // SDL_Log("%f %f", audio[i], audio[i + 1]);
        double current_zero_exact = it;
        if (it_next - it >= -it) {
            current_zero_exact = it_index + (-it / (it_next - it));
        } else continue;

        if (prev_zero != 0) {
            const int interval = it_index - prev_zero;
            const double interval_exact = current_zero_exact - prev_zero_exact;

            // отфильтровать интервал, проверить периодичность
            bool bad_zero_cross = false;
            {
                using std::fabs;
                float amplitude = 0.0f;
                FOR(range_idx, 0, interval) {
                    if (it_index + range_idx >= audio.size) break;
                    float val = fabs(audio[it_index + range_idx]);
                    if (val > amplitude) amplitude = val;
                }
                // TODO сделать чтобы примерно совпадаюшие, но разные по высоте считались хорошим интервалом, например предварительно нормализовать звук по амплитуде
                // TODO фильтровать совсем плохие данные, слишком тихие и тд
                const float tolerance = 0.3f * amplitude;

                double step = interval / interval_exact;
                FOR(range_idx, 0, interval) {
                    double x_left  = prev_zero_exact + range_idx * step;
                    double x_right = x_left + (current_zero_exact - prev_zero_exact);
                    if (int(x_right) + 1 >= audio.size) break;

                    double y_left  = get_spline_value(audio, x_left);
                    double y_right = get_spline_value(audio, x_right);
                    if (fabs(y_left - y_right) > tolerance) {
                        bad_zero_cross = true;
                        break;
                    }
                }
            }
            if (bad_zero_cross) continue;

            ++intervals_appearences[interval];
            intervals_sum_diff[interval] += (current_zero_exact - prev_zero_exact);
        }

        draw_vertical_line(context, 0x804040FF, (float) it_index / audio.size);
        prev_zero = it_index;
        prev_zero_exact = current_zero_exact;
    }

    int most_often_interval = 0;
    for (int i = 0; i < DEFAULT_SAMPLES; ++i) {
        if (intervals_appearences[i] > intervals_appearences[most_often_interval]) most_often_interval = i;
    }
    double most_often_interval_exact = intervals_sum_diff[most_often_interval] / intervals_appearences[most_often_interval]; // take avg

    double storage_intervals[DEFAULT_SAMPLES];
    Array <double> intervals(storage_intervals, 0, DEFAULT_SAMPLES);

    int storage_intervals_ids[DEFAULT_SAMPLES];
    Array <int> intervals_ids(storage_intervals_ids, 0, DEFAULT_SAMPLES);

    for (int i = 0; i < DEFAULT_SAMPLES; ++i) {
        if (intervals_appearences[i] > 0) {
            intervals.add(intervals_sum_diff[i] / intervals_appearences[i]);
            intervals_ids.add(i);
        }
    }

    static char buffer[100];
    int last_y = 0;
    if (print_intervals) {
        For (intervals) {
            SDL_Point line[2];
            line[0].x = (int) (it / DEFAULT_SAMPLES * context->win_w);
            line[0].y = 0;
            line[1].x = line[0].x;
            line[1].y = context->win_h;

            // TODO(den): рисовать тут горизонтальные линии красивые
            // SDL_SetRenderDrawColor(context->renderer, 0, 255, 0, 255); SDL_RenderDrawLines(context->renderer, line, 2);
            double Hz = context->spec.freq / it;
            SDL_snprintf(buffer, 100, "%d) %1.2f     %d", intervals_ids[it_index], Hz, intervals_appearences[intervals_ids[it_index]]);
            draw_text(context, buffer, line[0].x, last_y);

            last_y += 30;
        }
    }

    // SDL_Log("%d", most_often_interval);
    // double Hz = 1.0 / ((double) most_often_interval / context->spec.freq);
    double Hz = context->spec.freq / most_often_interval_exact;
    static char buffer_hz    [100];
    static char buffer_units [100];
    static char buffer_note  [100];
    if (buffer_hz[0] == 0 || Hz < 2000) {
        double units = hz_to_units(Hz);
        Note note = units_to_note(units);

        SDL_snprintf(buffer_hz,    100, "hz=%1.2f",    Hz);
        SDL_snprintf(buffer_units, 100, "units=%1.2f", units);
        SDL_snprintf(buffer_note,  100, "note=%s%d",   note.letter, note.octave);
        static double prev_Hz = 1.0;
        if (Hz > 0.1 && (fabs(Hz/prev_Hz) <= 0.99 || fabs(prev_Hz/Hz) <= 0.99)) {
            SDL_Log("%1.2f", Hz);
            prev_Hz = Hz;
        }
    }
    draw_text(context, buffer_hz,   context->win_w - 300, context->win_h - 180);
    draw_text(context, buffer_units,context->win_w - 300, context->win_h - 120);
    draw_text(context, buffer_note, context->win_w - 300, context->win_h - 60);
}

// ###########################################################################################
// ########################################### OLD ###########################################
// ###########################################################################################

void FFT(float * io_real, float * io_img, int size, bool forward);

struct Sorted_Peaks {
    int       *ids;
    float     *vals;
    int       n = 0;
    const int n_total;

    bool add_peak(int id, float value) {
        assert(ids);
        assert(vals);
        assert(n >= 0);
        assert(n <= n_total);
        assert(n_total >= 0);
        if (n_total == 0) return false;

        if (n == 0) {
            ids[0]  = id;
            vals[0] = value;
            n = 1;
            return true;
        }

        int insert_pos = 0;
        while (insert_pos < n && vals[insert_pos] > value) ++insert_pos;
        if (insert_pos == n_total) return false;

        if (n < n_total) {
            ids[n]  = ids[n-1];
            vals[n] = vals[n-1];
            ++n;
        }

        for (int i = n-1; i > insert_pos; --i) {
            ids[i]  = ids[i-1];
            vals[i] = vals[i-1];
        }

        ids[insert_pos]  = id;
        vals[insert_pos] = value;

        return true;
    }

    Sorted_Peaks(int *ids_storage, float *vals_storage, int storage_size) : ids(ids_storage), vals(vals_storage), n(0), n_total(storage_size) {}
};

void find_frequency_by_fft(Context *context, Array<float> const& audio) {
    static SDL_Point points      [DEFAULT_SAMPLES]; // size = width просто я хочу сделать resize поэтому не выношу width в константы
    static SDL_Point points_img  [DEFAULT_SAMPLES]; // size = width просто я хочу сделать resize поэтому не выношу width в константы
    static float storage_fft_result        [DEFAULT_SAMPLES];
    static float storage_fft_result_img    [DEFAULT_SAMPLES];
    static float storage_fft_result_module [DEFAULT_SAMPLES];
    // static float samples_orig              [DEFAULT_SAMPLES];
    float const* samples_orig = audio.data;

    for (int i = 0; i < DEFAULT_SAMPLES; ++i) {
        // samples_orig[i] = Get16bitAudioSample(&record->record_buffer[print_pos], context->spec.format); // TODO перенести spec в запись
        storage_fft_result[i] = samples_orig[i];
        storage_fft_result_img[i] = 0.0f;
        // print_pos += 2;
        // if (print_pos >= RECORD_BUFFER_SIZE) print_pos = 0;
    }
    // TODO Hamming window чтобы сгладить разрывы, если будет неочень результат
    FFT(storage_fft_result, storage_fft_result_img, DEFAULT_SAMPLES, true);
    for (int i = 0; i < DEFAULT_SAMPLES; ++i) {
        storage_fft_result_module[i] = sqrt(storage_fft_result[i] * storage_fft_result[i] + storage_fft_result_img[i] * storage_fft_result_img[i]);
    }

    const int fft_part_size = DEFAULT_SAMPLES / 2;
    float * fft_part     = storage_fft_result;
    float * fft_part_img = storage_fft_result_img;

    // TODO использовать Protrusion Value
    #define N_PEAKS 6
    int   peaks_ids[N_PEAKS];
    float peaks_vals[N_PEAKS];
    Sorted_Peaks peaks(peaks_ids, peaks_vals, N_PEAKS);

    calculate_peaks(&peaks, fft_part, fft_part_img, fft_part_size);

    // ЭТО БЛЯТЬ ВООБЩЕ НЕ РАБОТАЕТ НУЖНО ОТСЕИВАТЬ ПЛОХИЕ ПИКИ СОВСЕМ
    // define N_STORE_PEAKS 100
    // static float last_peaks[N_STORE_PEAKS];
    // static int   idx_to_write_new_peak = 0;
    // static int   peaks_remembered = 0;

    float actual_peak = 0.0f;
    { // https://www.codeproject.com/Articles/32172/FFT-Guitar-Tuner
        const double MIN_FREQ = 60;
        const double MAX_FREQ = 1300;

        double min_peak_value = 100000000000.0;
        int min_peak_index = 0;
        int min_optimal_interval = 0;

        const int verify_offset = 0;
                int verify_length = (int)(context->spec.freq / MIN_FREQ);
        for (int i = 0; i < peaks.n; i++) {
            int index = peaks.ids[i];
            if (index == 0) continue;
            int bin_start = fft_part_size / (index + 1);
            int bin_end   = fft_part_size / index;
            double peak_value = 1000000.0;
            int    interval   = 0;
            // scan bins frequencies/intervals
            {
                // distance between min and max range value can be big
                // limiting it to the fixed value
                const int max_steps = 30;
                int steps = bin_end - bin_start;
                if (steps > max_steps) steps = max_steps;
                else if (steps <= 0)   steps = 1;

                // trying all intervals in the range to find one with
                // smaller difference in signal waves
                for (int idx = 0; idx < steps; idx++)
                {
                    int tmp_interval = bin_start + (bin_end - bin_start) * idx / steps;

                    double sum = 0;
                    for (int j = 0; j < verify_length; j++)
                    {
                        //SDL_Log("taking samples_orig[%d] and samples_orig[%d]", index + j, index + j + interval);
                        double diff = samples_orig[index + j] - samples_orig[index + j + interval];
                        sum += diff * diff;
                    }
                    if (peak_value > sum)
                    {
                        peak_value = sum;
                        interval   = tmp_interval;
                    }
                }
            }

            if (peak_value < min_peak_value)
            {
                min_peak_value       = peak_value;
                min_peak_index       = index;
                min_optimal_interval = interval;
            }
        }
        if (min_optimal_interval > 0.00000001) {
            actual_peak = (float) context->spec.freq / min_optimal_interval;
        }
    }
#if 0
    if (peaks.n >= 3) {
        int i[3];
        i[0] = peaks.ids[0];
        i[1] = peaks.ids[1];
        i[2] = peaks.ids[2];

        std::sort(i, i+3);

        int i0 = i[0];
        int i1 = i[1];
        int i2 = i[2];

        float A = fft_part[i[0]];
        float B = fft_part[i[1]];
        float C = fft_part[i[2]];

        float top = (i1*i1 - i0*i0) * (C - B) - (i2*i2 - i1*i1) * (B - A);
        float bot =       (i1 - i0) * (C - B) -       (i2 - i1) * (B - A);
        float h = 0.5 * top / bot;
        last_peaks[idx_to_write_new_peak] = h;
        ++idx_to_write_new_peak;
        if (idx_to_write_new_peak == N_STORE_PEAKS) {
            idx_to_write_new_peak = 0;
        }
        if (peaks_remembered != N_STORE_PEAKS) {
            ++peaks_remembered;
        } else {
            float avg = 0;
            for (int i = 0; i < N_STORE_PEAKS; ++i) {
                avg += last_peaks[i];
            }
            avg /= N_STORE_PEAKS;
            actual_peak = avg;
        }


        // float x = 0.5 * (A - C) / (A + C - 2*B);
        // float k = ids[0] + x;
        // actual_peak = k * record_spec.freq / record_spec.samples;
    } else {
        ;
    }
#endif


    // data_to_points(points,      fft_part,     fft_part_size, context->win_w, context->win_h, BOTTOM);
    data_to_points(points,      storage_fft_result_module, fft_part_size, context->win_w, context->win_h, BOTTOM);
    // data_to_points(points_img,  fft_part_img, fft_part_size, context->win_w, context->win_h, TOP);


    SDL_Point max_lines[2 * N_PEAKS];
    //if (print_slower % 100 == 0) {
    for (int max_it = 0; max_it < peaks.n; ++max_it) {
        max_lines[max_it * 2].x     = (int) ((float) peaks.ids[max_it] / fft_part_size * context->win_w);
        max_lines[max_it * 2 + 1].x = max_lines[max_it * 2].x;
    }
    // std::sort(max_lines, max_lines + 2 * N_PEAKS, [](const SDL_Point& lhs, const SDL_Point& rhs) { return lhs.x < rhs.x; });
    for (int max_it = 0; max_it < peaks.n; ++max_it) {
        max_lines[max_it * 2 + (1 + max_it) % 2].y = points_img[max_lines[max_it * 2].x].y + 10;
        max_lines[max_it * 2 +       max_it % 2].y = points[max_lines[max_it * 2].x].y - 10;
    }
    //}
    SDL_SetRenderDrawColor(context->renderer, 0, 255, 0, 255); SDL_RenderDrawLines(context->renderer, points,      context->win_w);
    // SDL_SetRenderDrawColor(context->renderer, 255, 0, 0, 255); SDL_RenderDrawLines(context->renderer, points_img,  context->win_w);
    for (int max_it = 0; max_it < peaks.n; ++max_it) {
        SDL_SetRenderDrawColor(context->renderer, 255, 255, 255, 255); SDL_RenderDrawLines(context->renderer, &max_lines[max_it * 2], 2);
    }

    static char buffer[100];
    for (int i = 0; i < peaks.n; ++i) {
        // int Hz = 1.0 / ((double) (DEFAULT_SAMPLES - peaks.ids[i]) / context->spec.freq);
        int Hz = peaks.ids[i];
        draw_text(context, SDL_itoa(Hz, buffer, 10), (int) ((float) peaks.ids[i] / fft_part_size * context->win_w), 0);
        SDL_snprintf(buffer, 100, "%1.2f", peaks.vals[i]);
        draw_text(context, buffer, (int) ((float) peaks.ids[i] / fft_part_size * context->win_w), 40);
    }
    // SDL_snprintf(buffer, 100, "%1.2f", actual_peak);
    // draw_text(&context, buffer, 0, 100);
}

// COPYPASTA
// COPYPASTA
// COPYPASTA
// COPYPASTA
// COPYPASTA
// COPYPASTA
// COPYPASTA
// COPYPASTA
// COPYPASTA
// COPYPASTA
// COPYPASTA
// COPYPASTA
// COPYPASTA
// COPYPASTA
void visualize_wav(Audio_Wav const* audio_wav) {
    Context context; // TODO реструктурировать чтобы не хранить spec в двух местах
    context.spec = audio_wav->spec;

    context.window  = SDL_CreateWindow("Audio visualizer", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 1024, 768, SDL_WINDOW_OPENGL);
    if (!context.window) {SDL_Log("Could not create window: %s\n", SDL_GetError()); return;} defer {SDL_DestroyWindow(context.window);};
    SDL_GetWindowSize(context.window, &context.win_w, &context.win_h);

    context.renderer = SDL_CreateRenderer(context.window, -1, SDL_RENDERER_ACCELERATED);
    if (!context.renderer) {SDL_Log("Could not create rendered: %s\n", SDL_GetError()); return;} defer {SDL_DestroyRenderer(context.renderer);};

    if (TTF_Init() == -1) {SDL_Log("Could not init TTF: %s\n", SDL_GetError()); return;} defer {TTF_Quit();};

    char const *font_path = DATA_PATH "Poppins-Regular.ttf";
    context.font = TTF_OpenFont(font_path, 36);
    if (!context.font) {SDL_Log("Can't open font %s: %s", font_path, SDL_GetError()); return;} defer {TTF_CloseFont(context.font);};

    u64 samples_begin = 0;
    double dt = 0;
    std::chrono::high_resolution_clock clock;
    auto prev_time = clock.now();

    const double sample_time = (double) context.spec.samples / context.spec.freq;
    auto last_sample_moment = clock.now();
    // TODO синхронизация с мирофоном
    for(;;) {
        // SDL_Log("sample_tome %f", sample_time);
        auto current_time = clock.now();
        dt = std::chrono::duration<double>(current_time - prev_time).count();
        prev_time = current_time;

        double time_since_last_sample = std::chrono::duration<double>(current_time - last_sample_moment).count();
        if (time_since_last_sample >= 0.1) {
            // SDL_Log("time_since_last_sample = %f", time_since_last_sample);
            last_sample_moment = current_time;
            samples_begin += DEFAULT_SAMPLES * 2;
            if (samples_begin + DEFAULT_SAMPLES * 2 >= audio_wav->wav.size) {
                samples_begin = 0;
            }
        }

        ++print_slower;

        SDL_SetRenderDrawColor(context.renderer, 0, 0, 0, 255);
        SDL_RenderClear(context.renderer);

        static float samples_orig[DEFAULT_SAMPLES];
        for (int i = 0; i < DEFAULT_SAMPLES; ++i) {
          samples_orig[i] = Get16bitAudioSample(&audio_wav->wav[samples_begin + 2 * i], audio_wav->spec.format); // TODO перенести spec в запись
        }
        Array<float> audio(samples_orig, DEFAULT_SAMPLES);

        static SDL_Point points_orig [DEFAULT_SAMPLES]; // size = width просто я хочу сделать resize поэтому не выношу width в константы
        data_to_points(points_orig, samples_orig, DEFAULT_SAMPLES, context.win_w, context.win_h, MIDDLE);
        SDL_SetRenderDrawColor(context.renderer, 0, 0, 255, 255); SDL_RenderDrawLines(context.renderer, points_orig, context.win_w);

        // find_frequency_by_fft(&context, audio);
        find_frequency_by_zcr(&context, audio);

        SDL_RenderPresent(context.renderer);
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                return;
            }
        }
    }
}

void calculate_peaks(Sorted_Peaks * peaks, float * data, float * data_img, const int n_data) {
    float avg = 0.0f;
    for (int i = 0; i < n_data; ++i) {
        avg += data[i] * data[i] + data_img[i] * data_img[i];
    }
    avg /= n_data;

    const float MIN_PEAK_VALUE = avg;
    int idx = 0;
    while (idx < n_data) {
        while (idx < n_data) {
            float arr_val = data[idx] * data[idx] + data_img[idx] * data_img[idx];
            if (arr_val >= MIN_PEAK_VALUE) break;

            ++idx;
        }

        int   max_idx = 0;
        float max_val = 0.0f;
        while (idx < n_data) {
            float arr_val = data[idx] * data[idx] + data_img[idx] * data_img[idx];
            if (arr_val < MIN_PEAK_VALUE) break;

            if (arr_val > max_val) {
                max_val = arr_val;
                max_idx = idx;
            }

            ++idx;
        }

        peaks->add_peak(max_idx, max_val);
        ++idx;
    }
    if (print_slower % 200 == 0) {
        SDL_Log("avg = %f", avg);
    }
}

void multiply_complex(float* result_x, float* result_y, float l_x, float l_y, float r_x, float r_y) {
    *result_x = l_x * r_x - l_y * r_y;
    *result_y = l_x * r_y + l_y * r_x;
}
// https://github.com/d1vanov/Simple-FFT/blob/a0cc843ff36d33ad09c08674b9503614742ad0b9/include/simple_fft/fft_impl.hpp#L190
void FFT(float * io_real, float * io_img, int size, bool forward) {
    struct Complex { float x; float y; };
    assert(io_real);
    assert(io_img);

    const double local_pi = forward ? -M_PI : M_PI;
    size_t next, match;
    float sine;
    float delta;
    Complex mult, factor, product;

    for (size_t i = 1; i < size; i <<= 1) {
        next = i << 1;
        delta = float(local_pi / i);
        sine = (float) sin(0.5 * delta);
        mult.x = -2.0f * sine * sine;
        mult.y = sin(delta);
        factor.x = 1.0;
        factor.y = 0.0;

        for (size_t j = 0; j < i; ++j) {
            for (size_t k = j; k < size; k += next) {
                match = k + i;
                // product = fft_result[match] * factor;
                multiply_complex(&product.x, &product.y, io_real[match], io_img[match], factor.x, factor.y);

                io_real[match] = io_real[k] - product.x;
                io_img[match]  = io_img[k] - product.y;
                io_real[k] += product.x;
                io_img[k]  += product.y;
            }

            // factor = mult * factor + factor;
            Complex tmp;
            multiply_complex(&tmp.x, &tmp.y, mult.x, mult.y, factor.x, factor.y);
            factor.x = tmp.x + factor.x;
            factor.y = tmp.y + factor.y;
        }
    }
}

void play_from_mic(Devices* devices) {
    SDL_Log("play_from_mic(...)");
    Record record;
    Player_User_Data player_user_data;
    player_user_data.record = &record;

    SDL_AudioSpec play_spec, record_spec;
    set_default_spec(&play_spec);
    play_spec.callback = (SDL_AudioCallback)play_audio;
    play_spec.userdata = &player_user_data;
    set_default_spec(&record_spec);
    record_spec.callback = (SDL_AudioCallback)record_audio;
    record_spec.userdata = &record; // TODO typed callbacks

    SDL_Log("Opening devices...");
    SDL_AudioSpec tmp;
    SDL_AudioDeviceID play_id = SDL_OpenAudioDevice(devices->play_device, SDL_FALSE, &play_spec, &tmp, 0 /*SDL_AUDIO_ALLOW_FORMAT_CHANGE*/);
    if (play_id == 0) {SDL_Log("Failed to open audio play device: %s", SDL_GetError()); return; }
    defer {SDL_CloseAudioDevice(play_id);};
    if (tmp.format != play_spec.format) { SDL_Log("We didn't get audio play format: have %x, but given %x", tmp.format, play_spec.format); return; }

    SDL_AudioDeviceID record_id = SDL_OpenAudioDevice(devices->record_device, SDL_TRUE, &record_spec, &tmp, 0);
    if (record_id == 0) {SDL_Log("Failed to open audio record device: %s", SDL_GetError()); return; }
    defer {SDL_CloseAudioDevice(record_id);};
    if (tmp.format != record_spec.format) { SDL_Log("We didn't get audio record format: have %x, but given %x", tmp.format, record_spec.format); return; }

    SDL_Log("Recording...");
    SDL_PauseAudioDevice(record_id, 0);
    SDL_Log("Playing...");
    SDL_PauseAudioDevice(play_id, 0);
    SDL_Delay(1000000);
}

void play_audio(Player_User_Data * player_user_data, u8 * stream, int chunk_len) {
    if (player_user_data->position + chunk_len > RECORD_BUFFER_SIZE) {
        player_user_data->position = 0;
        SDL_Log(".");
    }
    u64 const end_pos   = player_user_data->record->current_writen_pos.load();
    u64 const start_pos = player_user_data->position;
    if (end_pos >= start_pos && start_pos + chunk_len > end_pos) {
        SDL_Log("skipped start = %d \t end = %d \t end_pos = %d", start_pos, start_pos + chunk_len, end_pos);
        return;
    }
    //SDL_Log("play audio %d %d", (int)start_pos, (int)end_pos);
    memcpy(stream, &player_user_data->record->record_buffer[start_pos], chunk_len);
    SDL_Log("wrote \t%d\t-\t%d", start_pos, start_pos + chunk_len);
    player_user_data->position = start_pos + chunk_len;
}

bool play_wav_from_drive(char const* const wav_file, char const* const play_device) {
    Player_State audio_data;
    Audio_Wav audio_wav;
    open_wav(&audio_wav, wav_file);
    defer {SDL_FreeWAV(audio_wav.wav.data);};

    audio_data.play_position = audio_wav.wav.data;
    audio_data.rest_length   = audio_wav.wav.size;

    SDL_AudioSpec tmp;
    SDL_AudioDeviceID playDevice = SDL_OpenAudioDevice(play_device, SDL_FALSE, &audio_wav.spec, &tmp, 0/*SDL_AUDIO_ALLOW_FORMAT_CHANGE*/);
    if (playDevice == 0) {SDL_Log("Failed to open audio device: %s", SDL_GetError()); return false; }
    defer {SDL_CloseAudioDevice(playDevice);};

    if (tmp.format != audio_wav.spec.format) { SDL_Log("We didn't get wav_spec audio format have %x, but given %x", tmp.format, audio_wav.spec.format); return false; }

    SDL_Log("Playing %s...", wav_file);
    SDL_PauseAudioDevice(playDevice, 0);
    while (audio_data.rest_length > 0) { // TODO: rest_length is shared, so technically it is data race?
        SDL_Delay(100);
    }

    return true;
}









