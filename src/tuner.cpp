
// нужно определить свой, просто не хочется чтобы этот путь зависил от того, в какой папке ты запустил бинарник
#define DATA_PATH "C:/Users/22599/projects/_game1/data/"

// глобальные переменные
float global_y_scale = 1.0f;

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

#include <assert.h>
#include "array.h"

#include "SDL.h"
#include "SDL_ttf.h"
#include <atomic>
#include <chrono>
#include <string>
#include <cmath>

using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

using u16 = uint16_t;
using s16 = int16_t;

using u64atom = std::atomic_uint64_t;

static constexpr int DEFAULT_SAMPLES = 4096;
static constexpr int RECORD_BUFFER_SIZE = DEFAULT_SAMPLES * 128;

struct Devices {
    char const* play_device   = NULL;
    char const* record_device = NULL;
};

struct Context {
    SDL_Window   *window;
    int win_w, win_h;
    SDL_Renderer *renderer;
    TTF_Font     *font;
    SDL_AudioSpec spec;
};

void set_default_spec(SDL_AudioSpec *spec);
bool get_record_device_from_cmd_args(char const **record_device, int n_args, char const* const* args);
void visualize_from_mic(char const* record_device);

int main(int const n_args, char ** const args) {
    if(SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO) == -1) {SDL_Log("Could not initialize SDL: %s.\n", SDL_GetError()); return -1;}
    defer {SDL_Quit();};

    char const* record_device = NULL;
    get_record_device_from_cmd_args(&record_device, n_args, args);
    if (!record_device) return -1;

    visualize_from_mic(record_device);

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

void print_record_devices(int n_record_devices) {
    SDL_Log("list of record devices:");
    for (int i = 0; i < n_record_devices; ++i) {
        SDL_Log("%d: %s", i, SDL_GetAudioDeviceName(i, SDL_TRUE));
    }
}

bool get_record_device_from_cmd_args(char const **record_device, int n_args, char const* const* args) {
    if (n_args != 1 && n_args != 2) {SDL_Log("need 0 or 1 arguments: %s idx_record", args[0]); return false;}

    int record_device_idx = -1;
    int n_record_devices = SDL_GetNumAudioDevices(SDL_TRUE);
    if (n_args == 1) {
        record_device_idx = 0;
    }
    if (n_args == 3) {
        record_device_idx = SDL_atoi(args[1]);
        if ((record_device_idx == 0 && args[1][0] != '0') || record_device_idx >= n_record_devices) {SDL_Log("you should choose recording device number from 0 to %d, but you entered %s", n_record_devices - 1, args[1]); return false;}
    }
    if (record_device_idx == -1) {SDL_Log("Pass the correct number as a second argument but not '%s'", args[1]); print_record_devices(n_record_devices); return false;}

    assert(record_device_idx < n_record_devices);
    *record_device = SDL_GetAudioDeviceName(record_device_idx, SDL_TRUE);

    SDL_Log("You chose recording device %d. %s", record_device_idx, record_device);
    print_record_devices(n_record_devices);
    return true;
}

// visualize_from_mic dependencies
struct Record {
    u64atom current_writen_pos = 0;
    Array<u8> record_buffer;
};

struct Audio_Wav {
    Array<u8>     wav;
    SDL_AudioSpec spec;
    char const*   name = NULL;
};

bool init_sdl(Context * context) {
    bool success = false;
    context->window  = SDL_CreateWindow("Audio visualizer", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 1024, 768, SDL_WINDOW_OPENGL); if (!context->window) {SDL_Log("Could not create window: %s\n", SDL_GetError());     return false;}
    SDL_GetWindowSize(context->window, &context->win_w, &context->win_h);
    defer {if (!success) SDL_DestroyWindow(context->window);};

    context->renderer = SDL_CreateRenderer(context->window, -1, SDL_RENDERER_ACCELERATED);                                                   if (!context->renderer) {SDL_Log("Could not create rendered: %s\n", SDL_GetError()); return false;}
    defer {if (!success) SDL_DestroyRenderer(context->renderer);};
    SDL_SetRenderDrawBlendMode(context->renderer, SDL_BLENDMODE_BLEND);

    int ttf_init_result = TTF_Init();                                                                                                        if (ttf_init_result == -1) {SDL_Log("Could not init TTF: %s\n", SDL_GetError());     return false;}
    defer {if (!success) TTF_Quit();};

    char const *font_path = DATA_PATH "Inconsolata-Regular.ttf";
    context->font = TTF_OpenFont(font_path, 36);                                                                                            if (!context->font) {SDL_Log("Can't open font %s: %s", font_path, SDL_GetError());    return false;}
    defer {if (!success) TTF_CloseFont(context->font);};

    success = true;
    return success;
}

void close_sdl(Context *context) {
    TTF_CloseFont(context->font);
    TTF_Quit();
    SDL_DestroyRenderer(context->renderer);
    SDL_DestroyWindow(context->window);
}

void record_audio(Record *record, u8 *stream, int chunk_len) {
    u64 start_pos = record->current_writen_pos.load();
    if (start_pos + chunk_len > RECORD_BUFFER_SIZE) start_pos = 0;

    memcpy(&record->record_buffer[start_pos], stream, chunk_len);
    record->current_writen_pos.store(start_pos + chunk_len);
}

enum Graph_Pos {
    TOP,
    BOTTOM,
    MIDDLE
};

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

void draw_audio(Context * context, Array<float> const& audio) {
    static SDL_Point points_orig [DEFAULT_SAMPLES]; // size = width просто я хочу сделать resize поэтому не выношу width в константы
    float const* samples_orig = audio.data;
    data_to_points(points_orig, samples_orig, audio.size, context->win_w, context->win_h, MIDDLE);
    SDL_SetRenderDrawColor(context->renderer, 0, 0, 255, 255); SDL_RenderDrawLines(context->renderer, points_orig, context->win_w);
    SDL_SetRenderDrawColor(context->renderer, 0, 255, 0, 255); SDL_RenderDrawLine(context->renderer, 0, context->win_h/2, context->win_w, context->win_h/2);
}

double absf(double x) {
    if (x >= 0.0) return x;
    else          return -x;
}

double get_spline_value(Array<float> const& audio, double x) {
    assert(x >= 0.0f);
    int    i = (int) x; // floor
    double e = x - i;
    return audio[i] + (audio[i+1] - audio[i]) * e;
}

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

void draw_vertical_line(Context *context, Color c, double pos) {
    assert(pos >= 0);
    assert(pos <= 1);
    SDL_SetRenderDrawColor(context->renderer, c.r, c.g, c.b, c.a);
    int screen_pos = (int) (pos * context->win_w);
    SDL_RenderDrawLine(context->renderer, screen_pos, 0, screen_pos, context->win_h);
}

void draw_text(Context *context, char const* text, int x, int y) {
    SDL_Rect message_rect;
    message_rect.x = x + 2;
    message_rect.y = y + 2;
    TTF_SizeText(context->font, text, &message_rect.w, &message_rect.h);

    {
        SDL_Surface* surfaceMessage = TTF_RenderText_Solid(context->font, text, {128, 128, 128}); if (!surfaceMessage) {SDL_Log("error: %s", SDL_GetError()); return;}
        defer {SDL_FreeSurface(surfaceMessage);};

        SDL_Texture* texture = SDL_CreateTextureFromSurface(context->renderer, surfaceMessage);   if (!texture) {SDL_Log("SDL_CreateTextureFromSurface: %s", SDL_GetError()); return;}
        defer {SDL_DestroyTexture(texture);};

        SDL_RenderCopy(context->renderer, texture, NULL, &message_rect);
    }

    message_rect.x = x;
    message_rect.y = y;

    {
        SDL_Surface* surfaceMessage = TTF_RenderText_Solid(context->font, text, {255, 255, 255}); if (!surfaceMessage) {SDL_Log("error: %s", SDL_GetError()); return;}
        defer {SDL_FreeSurface(surfaceMessage);};

        SDL_Texture* texture = SDL_CreateTextureFromSurface(context->renderer, surfaceMessage);   if (!texture) {SDL_Log("SDL_CreateTextureFromSurface: %s", SDL_GetError()); return;}
        defer {SDL_DestroyTexture(texture);};

        SDL_RenderCopy(context->renderer, texture, NULL, &message_rect);
    }
}

double hz_to_units(double hz) {
    // using log2
    static constexpr double log_of_2   = 1.0;
    static constexpr double log_of_440 = 8.781359713524659938;
    return 12 * (std::log2(hz) - log_of_440) / log_of_2;
}

struct Note {
    char const* letter;
    int         octave;
    double      cents;
};

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
    ret.cents  = (units - closest_note) * 100;
    return ret;
}

void find_frequency_by_zcr(Context * context, Array<float> const& audio) {
    static int    intervals_appearences[DEFAULT_SAMPLES];
    static double intervals_sum_diff   [DEFAULT_SAMPLES];
    for (int i = 0; i < DEFAULT_SAMPLES; ++i) intervals_appearences[i] = 0;
    for (int i = 0; i < DEFAULT_SAMPLES; ++i) intervals_sum_diff   [i] = 0;

    int    prev_zero       = 0;
    double prev_zero_exact = 0.0;

    For(audio) {
        if (it_index + 1 == audio.size) break;

        auto it_next = audio[it_index + 1];
        if (!(it < 0.0f && it_next >= 0.0f)) continue;

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
                float amplitude = 0.0f;
                FOR(range_idx, 0, interval) {
                    if (it_index + range_idx >= audio.size) break;
                    float val = (float)absf(audio[it_index + range_idx]);
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
                    if (absf(y_left - y_right) > tolerance) {
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

    double Hz = context->spec.freq / most_often_interval_exact;
    static char buffer_hz    [100];
    static char buffer_units [100];
    static char buffer_note  [100];
    if (buffer_hz[0] == 0 || Hz < 2000) {
        double units = hz_to_units(Hz);
        Note note = units_to_note(units);

        if (note.octave >= 0 && note.octave <= 9) {
            SDL_snprintf(buffer_hz,    100, "hz=%1.2f",    Hz);
            SDL_snprintf(buffer_units, 100, "units=%1.2f", units);
            SDL_snprintf(buffer_note,  100, "note=%s%s%d%s%s%.1f", note.letter, (note.letter[1] == 0 ? " " : ""), note.octave, (note.cents > 0 ? "+" : "-"), (abs(note.cents) < 10 ? "0" : ""), abs(note.cents));
            static double prev_Hz = 1.0;
            if (Hz > 0.1 && (absf(Hz/prev_Hz) <= 0.99 || absf(prev_Hz/Hz) <= 0.99)) {
                SDL_Log("%1.2f", Hz);
                prev_Hz = Hz;
            }
        }
    }
    draw_text(context, buffer_hz,   context->win_w - 400, context->win_h - 180);
    draw_text(context, buffer_units,context->win_w - 400, context->win_h - 120);
    draw_text(context, buffer_note, context->win_w - 400, context->win_h - 60);
}

void visualize_from_mic(char const* record_device) {
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
    record_id = SDL_OpenAudioDevice(record_device, SDL_TRUE, &context->spec, NULL, 0); if (record_id == 0) {SDL_Log("Failed to open audio record device: %s", SDL_GetError()); return; }
    defer {SDL_CloseAudioDevice(record_id);};

    SDL_PauseAudioDevice(record_id, 0);

    double dt = 0;
    std::chrono::high_resolution_clock clock;
    auto prev_time = clock.now();

    const double sample_time = (double) context->spec.samples / context->spec.freq;
    auto last_sample_moment = clock.now();
    for(;;) {
        auto current_time = clock.now();
        dt = std::chrono::duration<double>(current_time - prev_time).count();
        prev_time = current_time;

        double time_since_last_sample = std::chrono::duration<double>(current_time - last_sample_moment).count();
        if (time_since_last_sample >= sample_time) last_sample_moment = current_time;

        static u8 storage_samples_raw [DEFAULT_SAMPLES * 2];
        Audio_Wav wav_raw;
        wav_raw.spec = context->spec; // TODO убрать spec из контекста
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

        assert(context->spec.format == AUDIO_S16);
        for (int i = 0; i < DEFAULT_SAMPLES; ++i) {
            u16 val = (u16)storage_samples_raw[2*i] | ((u16)storage_samples_raw[2*i + 1] << 8);
            audio[i] = ((s16)val)/32768.0f;
        }

        SDL_SetRenderDrawColor(context->renderer, 0, 0, 0, 255);
        SDL_RenderClear(context->renderer);

        draw_audio(context, audio);
        find_frequency_by_zcr(context, audio);

        SDL_RenderPresent(context->renderer);

        // что-то немного странно все ивенты сразу съедать в кадре, а то ведь безумное шевеление мышкой может привести к зависанию?
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                return;
            } else if (event.type == SDL_KEYUP) {
                for (;;) {
                    if (SDL_PollEvent(&event)) {
                        // TODO сделать версию с keydown
                        if (event.type == SDL_KEYUP) {
                            break;
                        } else if (event.type == SDL_QUIT) {
                            return;
                        }
                    }
                }
            }
        }
    }
}
