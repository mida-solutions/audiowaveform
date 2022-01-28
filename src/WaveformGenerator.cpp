//------------------------------------------------------------------------------
//
// Copyright 2013-2021 BBC Research and Development
//
// Author: Chris Needham
//
// This file is part of Audio Waveform Image Generator.
//
// Audio Waveform Image Generator is free software: you can redistribute it
// and/or modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// Audio Waveform Image Generator is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// Audio Waveform Image Generator.  If not, see <http://www.gnu.org/licenses/>.
//
//------------------------------------------------------------------------------

#include "WaveformGenerator.h"
#include "Error.h"
#include "Log.h"
#include "WaveformBuffer.h"

#include <boost/format.hpp>

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <math.h> 

//------------------------------------------------------------------------------

ScaleFactor::~ScaleFactor()
{
}

//------------------------------------------------------------------------------

SamplesPerPixelScaleFactor::SamplesPerPixelScaleFactor(int samples_per_pixel) :
    samples_per_pixel_(samples_per_pixel)
{
}

//------------------------------------------------------------------------------

int SamplesPerPixelScaleFactor::getSamplesPerPixel(int /* sample_rate */) const
{
    return samples_per_pixel_;
}

//------------------------------------------------------------------------------

// Calculates samples_per_pixel such that the time range start_time to end_time
// fits the specified image width.

DurationScaleFactor::DurationScaleFactor(
    double start_time,
    double end_time,
    int width_pixels) :
    start_time_(start_time),
    end_time_(end_time),
    width_pixels_(width_pixels)
{
    if (end_time < start_time) {
        throwError("Invalid end time, must be greater than %1%", start_time);
    }

    if (width_pixels < 1) {
        throwError("Invalid image width: minimum 1");
    }
}

//------------------------------------------------------------------------------

int DurationScaleFactor::getSamplesPerPixel(const int sample_rate) const
{
    const double seconds = end_time_ - start_time_;

    const long long width_samples = static_cast<long long>(seconds * sample_rate);

    const int samples_per_pixel = static_cast<int>(width_samples / width_pixels_);

    return samples_per_pixel;
}

//------------------------------------------------------------------------------

PixelsPerSecondScaleFactor::PixelsPerSecondScaleFactor(int pixels_per_second) :
    pixels_per_second_(pixels_per_second)
{
    if (pixels_per_second_ <= 0) {
        throwError("Invalid pixels per second: must be greater than zero");
    }
}

//------------------------------------------------------------------------------

int PixelsPerSecondScaleFactor::getSamplesPerPixel(int sample_rate) const
{
    return sample_rate / pixels_per_second_;
}

//------------------------------------------------------------------------------

const int MAX_SAMPLE = std::numeric_limits<short>::max();
const int MIN_SAMPLE = std::numeric_limits<short>::min();

//------------------------------------------------------------------------------

WaveformGenerator::WaveformGenerator(
    WaveformBuffer& buffer,
    bool split_channels,
    const ScaleFactor& scale_factor) :
    buffer_(buffer),
    scale_factor_(scale_factor),
    split_channels_(split_channels),
    channels_(0),
    output_channels_(0),
    samples_per_pixel_(0),
    count_(0)
{
}

//------------------------------------------------------------------------------

bool WaveformGenerator::init(
    const int sample_rate,
    const int channels,
    const long /* frame_count */,
    const int /* buffer_size */)
{
    if (channels < 1 || channels > WaveformBuffer::MAX_CHANNELS) {
        log(Error) << "Cannot generate waveform data from audio file with "
                   << channels << " channels\n";
        return false;
    }

    channels_ = channels;

    samples_per_pixel_ = scale_factor_.getSamplesPerPixel(sample_rate);

    if (samples_per_pixel_ < 2) {
        log(Error) << "Invalid zoom: minimum 2\n";
        return false;
    }

    output_channels_ = split_channels_ ? channels : 1;

    buffer_.setSamplesPerPixel(samples_per_pixel_);
    buffer_.setSampleRate(sample_rate);
    buffer_.setChannels(output_channels_);

    log(Info) << "Generating waveform data...\n"
              << "Samples per pixel: " << samples_per_pixel_ << '\n'
              << "Input channels: " << channels_ << '\n'
              << "Output channels: " << output_channels_ << '\n';

    min_.resize(output_channels_, MAX_SAMPLE);
    max_.resize(output_channels_, MIN_SAMPLE);
    reset();

    return true;
}

//------------------------------------------------------------------------------

bool WaveformGenerator::shouldContinue() const
{
    return true;
}

//------------------------------------------------------------------------------

int WaveformGenerator::getSamplesPerPixel() const
{
    return samples_per_pixel_;
}

//------------------------------------------------------------------------------

void WaveformGenerator::reset()
{
    for (int channel = 0; channel < output_channels_; ++channel) {
        min_[channel] = MAX_SAMPLE;
        max_[channel] = MIN_SAMPLE;
    }

    count_ = 0;
}

//------------------------------------------------------------------------------

void WaveformGenerator::done()
{
    if (count_ > 0) {
        for (int channel = 0; channel < output_channels_; ++channel) {
            buffer_.appendSamples(
                static_cast<short>(min_[channel]),
                static_cast<short>(max_[channel])
            );
        }

        reset();
    }

    log(Info) << "Generated " << buffer_.getSize() << " points\n";
}

//------------------------------------------------------------------------------
double convert(short n, short bits) {
    double max = pow(2, bits - 1);
    double v = n < 0 ? n * max : n * (max - 1);
    double v2 = (max - 1) < v ? (max - 1) : v;
    return (-max) > v2 ? (-max) : v2;
}

// See BlockFile::CalcSummary in Audacity
//BEM qui calcolo
bool WaveformGenerator::process(
    const short* input_buffer,
    const int input_frame_count)
{

    short numPeaks = input_frame_count / samples_per_pixel_;
    int peaks[numPeaks*2];
    
    log(Info) << "BEM input_buffer:" << sizeof(input_buffer)<<"\n";
    log(Info) << "BEM input_frame_count:" << input_frame_count << "\n";
    log(Info) << "BEM channels_:" << channels_ << "\n";
    log(Info) << "BEM samples_per_pixel_:" << samples_per_pixel_ << "\n";




    //BEM TODO sistema bits
    short bits = 16;

    for (short i = 0; i < numPeaks; ++i) {
        short start = i * samples_per_pixel_;
        short end = (i + 1) * samples_per_pixel_ > input_frame_count ? input_frame_count : (i + 1) * samples_per_pixel_;


       /* for (short i = start; i < end; i++)
            log(Info) << input_buffer[i] << ",";
        log(Info) << "\n" << "\n" << "\n" << "\n" << "\n" << "\n";
        */

        short extrema_min = *(std::min_element(input_buffer + start, input_buffer + end));
        short extrema_max = *(std::max_element(input_buffer + start, input_buffer + end));
        //log(Info) << "extrema_min:"<< extrema_min<<"\n" << "extrema_max"<< extrema_max<<"\n" << "\n";
        int min = convert(extrema_min, bits);
        int max = convert(extrema_max, bits);
        peaks[i * 2] = extrema_min;
        peaks[i * 2 + 1] = extrema_max;
    }
    /*for (int i = 0; i < sizeof(peaks); i++)
        log(Info) << peaks[i]<< ",";
        */

    for (int i = 0; i< numPeaks; i++) {
        buffer_.appendSamples(
            peaks[i*2],peaks[i*2+1]
        );
    }

    //reset();
     
    return true;

    /*
    
    for (int i = 0; i < input_frame_count; ++i) {
        const int index = i * channels_;
        

        if (output_channels_ == 1) {
            // Sum samples from each input channel to make a single (mono) waveform
            int sample = 0;

            for (int channel = 0; channel < channels_; ++channel) {
                sample += input_buffer[index + channel];
            }

            sample /= channels_;

            // Avoid numeric overflow when converting to short
            if (sample > MAX_SAMPLE) {
                sample = MAX_SAMPLE;
            }
            else if (sample < MIN_SAMPLE) {
                sample = MIN_SAMPLE;
            }

            if (sample < min_[0]) {
                min_[0] = sample;
            }

            if (sample > max_[0]) {
                max_[0] = sample;
            }

            int bits = 16;
            if (extrema_min > sample)
                extrema_min = sample;
            if (extrema_max < sample)
                extrema_max = sample;

            min_[0] = convert(extrema_min, bits);
            max_[0] = convert(extrema_max, bits);
        }
        else {
            for (int channel = 0; channel < channels_; ++channel) {
                int sample = input_buffer[index + channel];
                if (extrema_min > sample)
                    extrema_min = sample;
                if (extrema_max < sample)
                    extrema_max = sample;
            }
            for (int channel = 0; channel < channels_; ++channel) {
                int sample = input_buffer[index + channel];

                // Avoid numeric overflow when converting to short
                if (sample > MAX_SAMPLE) {
                    sample = MAX_SAMPLE;
                }
                else if (sample < MIN_SAMPLE) {
                    sample = MIN_SAMPLE;
                }

                if (sample < min_[channel]) {
                    min_[channel] = sample;
                }

                if (sample > max_[channel]) {
                    max_[channel] = sample;
                }

                int bits = 16;

                min_[channel] = convert(extrema_min, bits);
                max_[channel] = convert(extrema_max, bits);
            }
        }

        if (++count_ == samples_per_pixel_) {
            for (int channel = 0; channel < output_channels_; ++channel) {
                buffer_.appendSamples(
                    static_cast<short>(min_[channel]),
                    static_cast<short>(max_[channel])
                );
            }

            reset();
        }
    }

    return true;*/
}

//------------------------------------------------------------------------------
