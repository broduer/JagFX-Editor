package app;

import javax.sound.sampled.*;
import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.border.LineBorder;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.plaf.basic.BasicSliderUI;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.geom.Path2D;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

class Buffer {
    private ByteArrayOutputStream bos;
    public int offset;
    public byte[] internalArray;

    public Buffer() {
        this.bos = new ByteArrayOutputStream();
        this.offset = 0;
    }

    public Buffer(byte[] data) {
        this.internalArray = data;
        this.offset = 0;
    }

    public byte[] toByteArray() {
        if (bos != null) {
            return bos.toByteArray();
        }
        return internalArray;
    }

    public void writeUnsignedByte(int value) {
        bos.write(value & 0xFF);
        offset++;
    }

    public void writeUnsignedShort(int value) {
        bos.write((value >>> 8) & 0xFF);
        bos.write(value & 0xFF);
        offset += 2;
    }

    public void writeInt(int value) {
        bos.write((byte)(value >> 24));
        bos.write((byte)(value >> 16));
        bos.write((byte)(value >> 8));
        bos.write((byte)value);
        offset += 4;
    }

    public void writeUShortSmart(int value) {
        if (value < 128) {
            writeUnsignedByte(value);
        } else {
            writeUnsignedByte((value >>> 8) | 0x80);
            writeUnsignedByte(value & 0xFF);
        }
    }

    public void writeShortSmart(int value) {
        if (value >= -64 && value < 64) {
            writeUnsignedByte(value + 64);
        } else {
            int unsignedShortValue = value + 49152;
            if (unsignedShortValue < 0 || unsignedShortValue >= 32768) {
                throw new IllegalArgumentException("Value out of range for ShortSmart: " + value + " (calculated unsigned: " + unsignedShortValue + ")");
            }
            writeUnsignedByte((unsignedShortValue >>> 8) | 0x80);
            writeUnsignedByte(unsignedShortValue & 0xFF);
        }
    }
    public int readUnsignedByte() {
        if (offset >= internalArray.length) {
            throw new ArrayIndexOutOfBoundsException("Attempt to read beyond buffer limit at offset " + offset + ". Buffer length: " + internalArray.length);
        }
        return internalArray[offset++] & 0xFF;
    }

    public int readUnsignedShort() {
        int b1 = readUnsignedByte();
        int b2 = readUnsignedByte();
        return (b1 << 8) | b2;
    }

    public int readInt() {
        int b1 = readUnsignedByte();
        int b2 = readUnsignedByte();
        int b3 = readUnsignedByte();
        int b4 = readUnsignedByte();
        return (b1 << 24) | (b2 << 16) | (b3 << 8) | b4;
    }

    public int readUShortSmart() {
        int val = readUnsignedByte();
        if (val < 128) {
            return val;
        } else {
            return ((val & 0x7F) << 8) | readUnsignedByte();
        }
    }

    public int readShortSmart() {
        int val = readUnsignedByte();
        if (val < 128) {
            return val - 64;
        } else {
            return ((val & 0x7F) << 8 | readUnsignedByte()) - 49152;
        }
    }
}

class ByteBufferUtils {

    public static void clearIntArray(int[] array, int offset, int length) {
        for (int i = offset; i < offset + length; i++) {
            array[i] = 0;
        }
    }

}

class SoundEnvelope {

    int segments;
    int[] durations;
    int[] phases;
    int start;
    int end;
    int form;
    int ticks;
    int phaseIndex;
    int step;
    int amplitude;
    int max;

    SoundEnvelope() {
        this.segments = 2;
        this.durations = new int[2];
        this.phases = new int[2];
        this.durations[1] = 65535;
        this.phases[1] = 65535;
        this.form = 0; // Default form to 0 (off)
    }

    final void decode(Buffer buffer) {
        this.form = buffer.readUnsignedByte();
        this.start = buffer.readInt();
        this.end = buffer.readInt();
        this.decodeSegments(buffer);
    }

    final void decodeSegments(Buffer buffer) {
        this.segments = buffer.readUnsignedByte();
        this.durations = new int[this.segments];
        this.phases = new int[this.segments];

        for (int segment = 0; segment < this.segments; ++segment) {
            this.durations[segment] = buffer.readUnsignedShort();
            this.phases[segment] = buffer.readUnsignedShort();
        }

    }

    final void reset() {
        this.ticks = 0;
        this.phaseIndex = 0;
        this.step = 0;
        this.amplitude = 0;
        this.max = 0;
    }

    final int doStep(int period) {
        if (this.max >= this.ticks) {
            this.amplitude = this.phases[this.phaseIndex++] << 15;
            if (this.phaseIndex >= this.segments) {
                this.phaseIndex = this.segments - 1;
            }

            this.ticks = (int)((double)this.durations[this.phaseIndex] / 65536.0D * (double) period);
            if (this.ticks > this.max) {
                this.step = ((this.phases[this.phaseIndex] << 15) - this.amplitude) / (this.ticks - this.max);
            }
        }

        this.amplitude += this.step;
        ++this.max;
        return this.amplitude - this.step >> 15;
    }

    public void encode(Buffer buffer) {
        buffer.writeUnsignedByte(this.form);
        buffer.writeInt(this.start);
        buffer.writeInt(this.end);
        this.encodeSegments(buffer);
    }

    final void encodeSegments(Buffer buffer) {
        buffer.writeUnsignedByte(this.segments);
        for (int segment = 0; segment < this.segments; ++segment) {
            buffer.writeUnsignedShort(this.durations[segment]);
            buffer.writeUnsignedShort(this.phases[segment]);
        }
    }
}

class SoundFilter {

    static float[][] minimizedCoefficients;
    static int[][] coefficients;
    static float forwardMinimizedCoefficientMultiplier;
    static int forwardMultiplier;
    int[] pairs;
    int[][][] phases;
    int[][][] magnitudes;
    int[] unity;

    static {
        minimizedCoefficients = new float[2][8];
        coefficients = new int[2][8];
    }

    SoundFilter() {
        this.pairs = new int[2];
        this.phases = new int[2][2][4];
        this.magnitudes = new int[2][2][4];
        this.unity = new int[2];
    }

    float adaptMagnitude(int direction, int i, float f) {
        float alpha = (float) this.magnitudes[direction][0][i] + f * (float)(this.magnitudes[direction][1][i] - this.magnitudes[direction][0][i]);
        alpha *= 0.0015258789F;
        return 1.0F - (float) Math.pow(10.0D, -alpha / 20.0F);
    }

    float adaptPhase(int direction, int i, float f) {
        float alpha = (float) this.phases[direction][0][i] + f * (float)(this.phases[direction][1][i] - this.phases[direction][0][i]);
        alpha *= 1.2207031E-4F;
        return normalize(alpha);
    }

    int compute(int direction, float f) {
        float magnitude;
        if (direction == 0) {
            magnitude = (float)this.unity[0] + (float)(this.unity[1] - this.unity[0]) * f;
            magnitude *= 0.0030517578F;
            forwardMinimizedCoefficientMultiplier = (float)Math.pow(0.1D, magnitude / 20.0F);
            forwardMultiplier = (int)(forwardMinimizedCoefficientMultiplier * 65536.0F);
        }

        if (this.pairs[direction] == 0) {
            return 0;
        } else {
            magnitude = this.adaptMagnitude(direction, 0, f);
            minimizedCoefficients[direction][0] = -2.0F * magnitude * (float) Math.cos(this.adaptPhase(direction, 0, f));
            minimizedCoefficients[direction][1] = magnitude * magnitude;

            float[] coefficientFloatArray;
            int pair;
            for (pair = 1; pair < this.pairs[direction]; ++pair) {
                magnitude = this.adaptMagnitude(direction, pair, f);
                float phase = -2.0F * magnitude * (float) Math.cos(this.adaptPhase(direction, pair, f));
                float coefficient = magnitude * magnitude;
                minimizedCoefficients[direction][pair * 2 + 1] = minimizedCoefficients[direction][pair * 2 - 1] * coefficient;
                minimizedCoefficients[direction][pair * 2] = minimizedCoefficients[direction][pair * 2 - 1] * phase + minimizedCoefficients[direction][pair * 2 - 2] * coefficient;

                for (int pair2 = pair * 2 - 1; pair2 >= 2; --pair2) {
                    coefficientFloatArray = minimizedCoefficients[direction];
                    coefficientFloatArray[pair2] += minimizedCoefficients[direction][pair2 - 1] * phase + minimizedCoefficients[direction][pair2 - 2] * coefficient;
                }

                coefficientFloatArray = minimizedCoefficients[direction];
                coefficientFloatArray[1] += minimizedCoefficients[direction][0] * phase + coefficient;
                coefficientFloatArray = minimizedCoefficients[direction];
                coefficientFloatArray[0] += phase;
            }

            if (direction == 0) {
                for (pair = 0; pair < this.pairs[0] * 2; ++pair) {
                    coefficientFloatArray = minimizedCoefficients[0];
                    coefficientFloatArray[pair] *= forwardMinimizedCoefficientMultiplier;
                }
            }

            for (pair = 0; pair < this.pairs[direction] * 2; ++pair) {
                coefficients[direction][pair] = (int)(minimizedCoefficients[direction][pair] * 65536.0F);
            }

            return this.pairs[direction] * 2;
        }
    }

    final void decode(Buffer buffer, SoundEnvelope envelope) {
        int count = buffer.readUnsignedByte();
        this.pairs[0] = count >> 4;
        this.pairs[1] = count & 15;

        if (count != 0) {
            this.unity[0] = buffer.readUnsignedShort();
            this.unity[1] = buffer.readUnsignedShort();
            int migrated = buffer.readUnsignedByte();

            int direction;
            int pair;
            for (direction = 0; direction < 2; ++direction) {
                for (pair = 0; pair < this.pairs[direction]; ++pair) {
                    this.phases[direction][0][pair] = buffer.readUnsignedShort();
                    this.magnitudes[direction][0][pair] = buffer.readUnsignedShort();
                }
            }

            for (direction = 0; direction < 2; ++direction) {
                for (pair = 0; pair < this.pairs[direction]; ++pair) {
                    if ((migrated & (1 << (direction * 4 + pair))) != 0) {
                        this.phases[direction][1][pair] = buffer.readUnsignedShort();
                        this.magnitudes[direction][1][pair] = buffer.readUnsignedShort();
                    } else {
                        this.phases[direction][1][pair] = this.phases[direction][0][pair];
                        this.magnitudes[direction][1][pair] = this.magnitudes[direction][0][pair];
                    }
                }
            }

            if (migrated != 0 || this.unity[1] != this.unity[0]) {
                envelope.decodeSegments(buffer);
            }
        } else {
            int[] unityArray = this.unity;
            this.unity[1] = 0;
            unityArray[0] = 0;
        }
    }

    static float normalize(float alpha) {
        float f = 32.703197F * (float) Math.pow(2.0D, alpha);
        return f * 3.1415927F / 11025.0F;
    }

    /**
     * Encodes the SoundFilter parameters into the provided Buffer.
     * This method is designed to be the inverse of the decode method.
     *
     * @param envelope The SoundEnvelope associated with this filter, used for conditional segment encoding.
     */
    public void encode(Buffer buffer, SoundEnvelope envelope) {
        // Clamp pairs values when writing to ensure compatibility with decoders
        // that expect max 4 pairs.
        int clampedPairs0 = Math.min(this.pairs[0], 4);
        int clampedPairs1 = Math.min(this.pairs[1], 4);
        int count = (clampedPairs0 << 4) | (clampedPairs1 & 0xF);
        buffer.writeUnsignedByte(count);

        if (count != 0) {
            buffer.writeUnsignedShort(this.unity[0]);
            buffer.writeUnsignedShort(this.unity[1]);

            int migrated = getMigrated(clampedPairs0, clampedPairs1);
            buffer.writeUnsignedByte(migrated);

            for (int direction = 0; direction < 2; ++direction) {
                for (int pair = 0; pair < (direction == 0 ? clampedPairs0 : clampedPairs1); ++pair) {
                    buffer.writeUnsignedShort(this.phases[direction][0][pair]);
                    buffer.writeUnsignedShort(this.magnitudes[direction][0][pair]);
                }
            }

            for (int direction = 0; direction < 2; ++direction) {
                for (int pair = 0; pair < (direction == 0 ? clampedPairs0 : clampedPairs1); ++pair) {
                    if ((migrated & (1 << (direction * 4 + pair))) != 0) {
                        buffer.writeUnsignedShort(this.phases[direction][1][pair]);
                        buffer.writeUnsignedShort(this.magnitudes[direction][1][pair]);
                    }
                }
            }

            if (migrated != 0 || this.unity[1] != this.unity[0]) {
                envelope.encodeSegments(buffer);
            }
        }
    }

    private int getMigrated(int clampedPairs0, int clampedPairs1) {
        int migrated = 0;
        for (int direction = 0; direction < 2; ++direction) {
            for (int pair = 0; pair < (direction == 0 ? clampedPairs0 : clampedPairs1); ++pair) {
                if (this.phases[direction][1][pair] != this.phases[direction][0][pair] ||
                        this.magnitudes[direction][1][pair] != this.magnitudes[direction][0][pair]) {
                    migrated |= (1 << (direction * 4 + pair));
                }
            }
        }
        return migrated;
    }
}

// Integrated SoundTone class
class SoundTone {

    static int[] toneSamples;
    static int[] toneNoise;
    static int[] toneSine;
    static int[] tonePhases;
    static int[] toneDelays;
    static int[] toneVolumeSteps;
    static int[] tonePitchSteps;
    static int[] tonePitchBaseSteps;
    SoundEnvelope pitch; // Always present
    SoundEnvelope volume; // Always present
    SoundEnvelope pitchModifier = null; // Optional
    SoundEnvelope pitchModifierAmplitude = null; // Optional
    SoundEnvelope volumeMultiplier = null; // Optional
    SoundEnvelope volumeMultiplierAmplitude = null; // Optional
    SoundEnvelope release = null; // Optional
    SoundEnvelope attack = null; // Optional
    int[] oscillatorVolume;
    int[] oscillatorPitch;
    int[] oscillatorDelays;
    int delayTime;
    int delayDecay;
    SoundFilter filter; // Always present
    SoundEnvelope filterEnvelope; // Always present
    public int duration;
    public int offset;

    static {
        toneNoise = new int[32768];
        Random randomID = new Random(0L);

        int toneID;
        for (toneID = 0; toneID < 32768; ++toneID) {
            toneNoise[toneID] = (randomID.nextInt() & 2) - 1;
        }

        toneSine = new int[32768];

        for (toneID = 0; toneID < 32768; ++toneID) {
            toneSine[toneID] = (int)(Math.sin((double) toneID / 5215.1903D) * 16384.0D);
        }

        toneSamples = new int[220500];
        tonePhases = new int[5];
        toneDelays = new int[5];
        toneVolumeSteps = new int[5];
        tonePitchSteps = new int[5];
        tonePitchBaseSteps = new int[5];
    }

    public SoundTone() {
        this.oscillatorVolume = new int[]{0, 0, 0, 0, 0};
        this.oscillatorPitch = new int[]{0, 0, 0, 0, 0};
        this.oscillatorDelays = new int[]{0, 0, 0, 0, 0};
        this.delayTime = 0;
        this.delayDecay = 100;
        this.duration = 500;
        this.offset = 0;

        // Initialize mandatory envelopes and filter
        this.pitch = new SoundEnvelope();
        this.volume = new SoundEnvelope();
        this.filter = new SoundFilter();
        this.filterEnvelope = new SoundEnvelope();
    }

    public final int[] synthesize(int steps, int tones) {
        ByteBufferUtils.clearIntArray(toneSamples, 0, steps);
        if (tones >= 10) { // This condition seems arbitrary, should probably be based on duration
            double duration = (double) steps / ((double) tones + 0.0D);
            this.pitch.reset();
            this.volume.reset();
            int pitchModulationStep = 0;
            int pitchModulationBaseStep = 0;
            int pitchModulationPhase = 0;
            if (this.pitchModifier != null) {
                this.pitchModifier.reset();
                this.pitchModifierAmplitude.reset();
                pitchModulationStep = (int) ((double) (this.pitchModifier.end - this.pitchModifier.start) * 32.768D / duration);
                pitchModulationBaseStep = (int) ((double) this.pitchModifier.start * 32.768D / duration);
            }

            int volumeModulationStep = 0;
            int volumeModulationBaseStep = 0;
            int volumeModulationPhase = 0;
            if (this.volumeMultiplier != null) {
                this.volumeMultiplier.reset();
                this.volumeMultiplierAmplitude.reset();
                volumeModulationStep = (int) ((double) (this.volumeMultiplier.end - this.volumeMultiplier.start) * 32.768D / duration);
                volumeModulationBaseStep = (int) ((double) this.volumeMultiplier.start * 32.768D / duration);
            }

            int step;
            for (step = 0; step < 5; ++step) {
                if (this.oscillatorVolume[step] != 0) {
                    tonePhases[step] = 0;
                    toneDelays[step] = (int) ((double) this.oscillatorDelays[step] * duration);
                    toneVolumeSteps[step] = (this.oscillatorVolume[step] << 14) / 100;
                    tonePitchSteps[step] = (int) ((double) (this.pitch.end - this.pitch.start) * 32.768D * Math.pow(1.0057929410678534D, this.oscillatorPitch[step]) / duration);
                    tonePitchBaseSteps[step] = (int) ((double) this.pitch.start * 32.768D / duration);
                }
            }

            int pitchChange;
            int volumeChange;
            int volumeMultiplierChange;
            int volumeMultiplierAmplitudeChange;
            int[] samples;
            for (step = 0; step < steps; ++step) {
                pitchChange = this.pitch.doStep(steps);
                volumeChange = this.volume.doStep(steps);
                if (this.pitchModifier != null) {
                    volumeMultiplierChange = this.pitchModifier.doStep(steps);
                    volumeMultiplierAmplitudeChange = this.pitchModifierAmplitude.doStep(steps);
                    pitchChange += this.evaluateWave(pitchModulationPhase, volumeMultiplierAmplitudeChange, this.pitchModifier.form) >> 1;
                    pitchModulationPhase = pitchModulationPhase + pitchModulationBaseStep + (volumeMultiplierChange * pitchModulationStep >> 16);
                }

                if (this.volumeMultiplier != null) {
                    volumeMultiplierChange = this.volumeMultiplier.doStep(steps);
                    volumeMultiplierAmplitudeChange = this.volumeMultiplierAmplitude.doStep(steps);
                    volumeChange = volumeChange * ((this.evaluateWave(volumeModulationPhase, volumeMultiplierAmplitudeChange, this.volumeMultiplier.form) >> 1) + 32768) >> 15;
                    volumeModulationPhase = volumeModulationPhase + volumeModulationBaseStep + (volumeMultiplierChange * volumeModulationStep >> 16);
                }

                for (volumeMultiplierChange = 0; volumeMultiplierChange < 5; ++volumeMultiplierChange) {
                    if (this.oscillatorVolume[volumeMultiplierChange] != 0) {
                        volumeMultiplierAmplitudeChange = toneDelays[volumeMultiplierChange] + step;
                        if (volumeMultiplierAmplitudeChange < steps) {
                            samples = toneSamples;
                            samples[volumeMultiplierAmplitudeChange] += this.evaluateWave(tonePhases[volumeMultiplierChange], volumeChange * toneVolumeSteps[volumeMultiplierChange] >> 15, this.pitch.form);
                            samples = tonePhases;
                            samples[volumeMultiplierChange] += (pitchChange * tonePitchSteps[volumeMultiplierChange] >> 16) + tonePitchBaseSteps[volumeMultiplierChange];
                        }
                    }
                }
            }

            int volumeAttackAmplitudeChange;
            if (this.release != null) {
                this.release.reset();
                this.attack.reset();
                step = 0;
                boolean muted = true;

                for (volumeMultiplierChange = 0; volumeMultiplierChange < steps; ++volumeMultiplierChange) {
                    volumeMultiplierAmplitudeChange = this.release.doStep(steps);
                    volumeAttackAmplitudeChange = this.attack.doStep(steps);
                    if (muted) {
                        pitchChange = (volumeMultiplierAmplitudeChange * (this.release.end - this.release.start) >> 8) + this.release.start;
                    } else {
                        pitchChange = (volumeAttackAmplitudeChange * (this.release.end - this.release.start) >> 8) + this.release.start;
                    }

                    step += 256;
                    if (step >= pitchChange) {
                        step = 0;
                        muted = !muted;
                    }

                    if (muted) {
                        toneSamples[volumeMultiplierChange] = 0;
                    }
                }
            }

            if (this.delayTime > 0 && this.delayDecay > 0) {
                step = (int) ((double) this.delayTime * duration);

                for (pitchChange = step; pitchChange < steps; ++pitchChange) {
                    samples = toneSamples;
                    samples[pitchChange] += toneSamples[pitchChange - step] * this.delayDecay / 100;
                }
            }

            if (this.filter.pairs[0] > 0 || this.filter.pairs[1] > 0) {
                this.filterEnvelope.reset();
                step = this.filterEnvelope.doStep(steps + 1);
                pitchChange = this.filter.compute(0, (float) step / 65536.0F);
                volumeChange = this.filter.compute(1, (float) step / 65536.0F);
                if (steps >= pitchChange + volumeChange) {
                    volumeMultiplierChange = 0;
                    volumeMultiplierAmplitudeChange = Math.min(volumeChange, steps - pitchChange);
                    int var17;
                    while (volumeMultiplierChange < volumeMultiplierAmplitudeChange) {
                        volumeAttackAmplitudeChange = (int) ((long) toneSamples[volumeMultiplierChange + pitchChange] * (long) SoundFilter.forwardMultiplier >> 16);

                        for (var17 = 0; var17 < pitchChange; ++var17) {
                            volumeAttackAmplitudeChange += (int) ((long) toneSamples[volumeMultiplierChange + pitchChange - 1 - var17] * (long) SoundFilter.coefficients[0][var17] >> 16);
                        }

                        for (var17 = 0; var17 < volumeMultiplierChange; ++var17) {
                            volumeAttackAmplitudeChange -= (int) ((long) toneSamples[volumeMultiplierChange - 1 - var17] * (long) SoundFilter.coefficients[1][var17] >> 16);
                        }

                        toneSamples[volumeMultiplierChange] = volumeAttackAmplitudeChange;
                        step = this.filterEnvelope.doStep(steps + 1);
                        ++volumeMultiplierChange;
                    }

                    volumeMultiplierAmplitudeChange = 128;

                    while (true) {
                        if (volumeMultiplierAmplitudeChange > steps - pitchChange) {
                            volumeMultiplierAmplitudeChange = steps - pitchChange;
                        }

                        int var18;
                        while (volumeMultiplierChange < volumeMultiplierAmplitudeChange) {
                            var17 = (int) ((long) toneSamples[volumeMultiplierChange + pitchChange] * (long) SoundFilter.forwardMultiplier >> 16);

                            for (var18 = 0; var18 < pitchChange; ++var18) {
                                var17 += (int) ((long) toneSamples[volumeMultiplierChange + pitchChange - 1 - var18] * (long) SoundFilter.coefficients[0][var18] >> 16);
                            }

                            for (var18 = 0; var18 < volumeChange; ++var18) {
                                var17 -= (int) ((long) toneSamples[volumeMultiplierChange - 1 - var18] * (long) SoundFilter.coefficients[1][var18] >> 16);
                            }

                            toneSamples[volumeMultiplierChange] = var17;
                            step = this.filterEnvelope.doStep(steps + 1);
                            ++volumeMultiplierChange;
                        }

                        if (volumeMultiplierChange >= steps - pitchChange) {
                            while (volumeMultiplierChange < steps) {
                                var17 = 0;

                                for (var18 = volumeMultiplierChange + pitchChange - steps; var18 < pitchChange; ++var18) {
                                    var17 += (int) ((long) toneSamples[volumeMultiplierChange + pitchChange - 1 - var18] * (long) SoundFilter.coefficients[0][var18] >> 16);
                                }

                                for (var18 = 0; var18 < volumeChange; ++var18) {
                                    var17 -= (int) ((long) toneSamples[volumeMultiplierChange - 1 - var18] * (long) SoundFilter.coefficients[1][var18] >> 16);
                                }

                                toneSamples[volumeMultiplierChange] = var17;
                                this.filterEnvelope.doStep(steps + 1);
                                ++volumeMultiplierChange;
                            }
                            break;
                        }

                        pitchChange = this.filter.compute(0, (float) step / 65536.0F);
                        volumeChange = this.filter.compute(1, (float) step / 65536.0F);
                        volumeMultiplierAmplitudeChange += 128;
                    }
                }
            }

            for (step = 0; step < steps; ++step) {
                if (toneSamples[step] < -32768) {
                    toneSamples[step] = -32768;
                }

                if (toneSamples[step] > 32767) {
                    toneSamples[step] = 32767;
                }
            }

        }
        return toneSamples;
    }

    final int evaluateWave(int var1, int var2, int var3) {
        if (var3 == 1) {
            return (var1 & 32767) < 16384 ? var2 : -var2;
        } else if (var3 == 2) {
            return toneSine[var1 & 32767] * var2 >> 14;
        } else if (var3 == 3) {
            return (var2 * (var1 & 32767) >> 14) - var2;
        } else {
            return var3 == 4 ? var2 * toneNoise[var1 / 2607 & 32767] : 0;
        }
    }

    public final void decode(Buffer var1) {
        this.pitch = new SoundEnvelope();
        this.pitch.decode(var1);
        this.volume = new SoundEnvelope();
        this.volume.decode(var1);

        int var2 = var1.readUnsignedByte();
        if (var2 != 0) {
            --var1.offset;
            this.pitchModifier = new SoundEnvelope();
            this.pitchModifier.decode(var1);
            this.pitchModifierAmplitude = new SoundEnvelope();
            this.pitchModifierAmplitude.decode(var1);
        } else {
            this.pitchModifier = null;
            this.pitchModifierAmplitude = null;
        }

        var2 = var1.readUnsignedByte();
        if (var2 != 0) {
            --var1.offset;
            this.volumeMultiplier = new SoundEnvelope();
            this.volumeMultiplier.decode(var1);
            this.volumeMultiplierAmplitude = new SoundEnvelope();
            this.volumeMultiplierAmplitude.decode(var1);
        } else {
            this.volumeMultiplier = null;
            this.volumeMultiplierAmplitude = null;
        }

        var2 = var1.readUnsignedByte();
        if (var2 != 0) {
            --var1.offset;
            this.release = new SoundEnvelope();
            this.release.decode(var1);
            this.attack = new SoundEnvelope();
            this.attack.decode(var1);
        } else {
            this.release = null;
            this.attack = null;
        }

        for (int var3 = 0; var3 < 5; ++var3) {
            int var4 = var1.readUShortSmart();
            if (var4 == 0) {
                for (int i = var3; i < 5; i++) {
                    this.oscillatorVolume[i] = 0;
                    this.oscillatorPitch[i] = 0;
                    this.oscillatorDelays[i] = 0;
                }
                break;
            }
            this.oscillatorVolume[var3] = var4;
            this.oscillatorPitch[var3] = var1.readShortSmart();
            this.oscillatorDelays[var3] = var1.readUShortSmart();
        }

        this.delayTime = var1.readUShortSmart();
        this.delayDecay = var1.readUShortSmart();
        this.duration = var1.readUnsignedShort();
        this.offset = var1.readUnsignedShort();
        this.filter = new SoundFilter(); // Always instantiate filter
        this.filterEnvelope = new SoundEnvelope(); // Always instantiate filterEnvelope
        this.filter.decode(var1, this.filterEnvelope);
    }

    public void encode(Buffer buffer) {
        this.pitch.encode(buffer); // Always present, never null
        this.volume.encode(buffer); // Always present, never null

        // Pitch Modifier and Pitch Modifier Amplitude
        // Only encode if the object exists AND its form is non-zero
        if (this.pitchModifier != null && this.pitchModifier.form != 0) {
            this.pitchModifier.encode(buffer);
            this.pitchModifierAmplitude.encode(buffer);
        } else {
            buffer.writeUnsignedByte(0); // Flag for absent pitchModifier (its form byte)
        }

        // Volume Multiplier and Volume Multiplier Amplitude
        if (this.volumeMultiplier != null && this.volumeMultiplier.form != 0) {
            this.volumeMultiplier.encode(buffer);
            this.volumeMultiplierAmplitude.encode(buffer);
        } else {
            buffer.writeUnsignedByte(0);
        }

        // Release and Attack
        if (this.release != null && this.release.form != 0) {
            this.release.encode(buffer);
            this.attack.encode(buffer);
        } else {
            buffer.writeUnsignedByte(0);
        }

        // Oscillator Data (now limited to 5 as per original SoundTone)
        boolean terminatorWritten = false; // Flag to indicate if a 0 terminator has been written
        for (int i = 0; i < 5; ++i) { // Loop only 5 times for the 5 oscillators
            if (this.oscillatorVolume[i] != 0 || this.oscillatorPitch[i] != 0 || this.oscillatorDelays[i] != 0) {
                buffer.writeUShortSmart(this.oscillatorVolume[i]); // Write volume
                buffer.writeShortSmart(this.oscillatorPitch[i]); // Write pitch
                buffer.writeUShortSmart(this.oscillatorDelays[i]); // Write delay
            } else {
                buffer.writeUShortSmart(0); // Terminator for the oscillator list
                terminatorWritten = true;
                break; // Stop writing further oscillators
            }
        }
        // If all 5 oscillators had non-zero values, we still need to write a 0 terminator
        if (!terminatorWritten) {
            buffer.writeUShortSmart(0);
        }


        buffer.writeUShortSmart(this.delayTime); // Changed to UShortSmart as per provided SoundTone.java
        buffer.writeUShortSmart(this.delayDecay); // Changed to UShortSmart as per provided SoundTone.java
        buffer.writeUnsignedShort(this.duration);
        buffer.writeUnsignedShort(this.offset);

        // Filter and Filter Envelope
        this.filter.encode(buffer, this.filterEnvelope);
    }
}

class AudioDataSource {
    byte[] audioDataByte; // Changed to byte[] as mix() now returns byte[]
    int sampleRate;
    int startSample; // Added to match constructor usage in SoundEffect
    int endSample;   // Added to match constructor usage in SoundEffect

    // Constructor for byte[]
    public AudioDataSource(int sampleRate, byte[] audioDataByte, int start, int end) {
        this.sampleRate = sampleRate;
        this.audioDataByte = audioDataByte;
        this.startSample = start;
        this.endSample = end;
    }
}

// Integrated SoundEffect class
class SoundEffect {

    SoundTone[] soundTones;
    int start;
    int end;
    byte[] rawAudioData;

    public SoundEffect(Buffer var1) {
        this.soundTones = new SoundTone[10];

        for (int var2 = 0; var2 < 10; ++var2) {
            int var3 = var1.readUnsignedByte();
            if (var3 != 0) {
                --var1.offset;
                this.soundTones[var2] = new SoundTone();
                this.soundTones[var2].decode(var1);
            }
        }

        this.start = var1.readUnsignedShort();
        this.end = var1.readUnsignedShort();
    }

    // Default constructor for creating a new SoundEffect
    public SoundEffect() {
        this.soundTones = new SoundTone[10];
        this.start = 0;
        this.end = 0;
    }

    final byte[] mix() {
        int var1 = 0;

        int var2;
        for (var2 = 0; var2 < 10; ++var2) {
            if (this.soundTones[var2] != null && this.soundTones[var2].duration + this.soundTones[var2].offset > var1) {
                var1 = this.soundTones[var2].duration + this.soundTones[var2].offset;
            }
        }

        if (var1 == 0) {
            return new byte[0];
        } else {
            var2 = var1 * 22050 / 1000;
            byte[] var3 = new byte[var2];

            for (int var4 = 0; var4 < 10; ++var4) {
                if (this.soundTones[var4] != null) {
                    int var5 = this.soundTones[var4].duration * 22050 / 1000;
                    int var6 = this.soundTones[var4].offset * 22050 / 1000;
                    int[] var7 = this.soundTones[var4].synthesize(var5, this.soundTones[var4].duration);

                    for (int var8 = 0; var8 < var5; ++var8) {
                        int var9 = (var7[var8] >> 8) + var3[var8 + var6];
                        if ((var9 + 128 & -256) != 0) {
                            var9 = var9 >> 31 ^ 127;
                        }

                        var3[var8 + var6] = (byte)var9;
                    }
                }
            }

            return var3;
        }
    }

    public AudioDataSource toRawSound() {
        new Thread(() -> rawAudioData = this.mix()).start();
        return new AudioDataSource(22050, rawAudioData, this.start * 22050 / 1000, this.end * 22050 / 1000);
    }

    public byte[] encode() {
        Buffer buffer = new Buffer();
        for (int i = 0; i < 10; ++i) {
            if (this.soundTones[i] != null) {
                this.soundTones[i].encode(buffer);
            } else {
                buffer.writeUnsignedByte(0);
            }
        }
        buffer.writeUnsignedShort(this.start);
        buffer.writeUnsignedShort(this.end);
        return buffer.toByteArray();
    }
}


public class SoundEffectEditor extends JFrame {

    private static final Color BG_DEFAULT_FRAME = new Color(30, 30, 30); // Dark background
    private static final Color PANEL_BG = new Color(45, 45, 45); // Slightly lighter dark for panels
    private static final Color BORDER_COLOR = new Color(60, 60, 60); // Dark grey for borders
    private static final Color TEXT_COLOR_LIGHT = new Color(220, 220, 220); // Light text for contrast
    private static final Color HIGHLIGHT_COLOR_GREEN = new Color(0, 255, 0); // Vibrant green for highlights
    private static final Color BUTTON_BG = new Color(70, 70, 70); // Dark button grey
    private static final Color BUTTON_BLUE = new Color(50, 150, 250); // Brighter blue for actions
    private static final Color BUTTON_RED = new Color(200, 50, 50); // Red for destructive actions
    private static final Color BUTTON_GREEN = new Color(50, 200, 50); // Green for positive actions
    private static final Color SLIDER_THUMB_COLOR = HIGHLIGHT_COLOR_GREEN; // Vibrant green thumb
    private static final Color SLIDER_TRACK_COLOR = new Color(90, 90, 90); // Darker track color
    private static final Color GRAPH_BG = Color.BLACK; // Black background for graphs
    private static final Color GRAPH_GRID = new Color(60, 60, 60); // Dark grey grid
    private static final Color GRAPH_LINE = HIGHLIGHT_COLOR_GREEN; // Vibrant green graph line
    private static final Color WAVEFORM_BG = Color.BLACK; // Black background for waveform
    private static final Color WAVEFORM_COLOR = HIGHLIGHT_COLOR_GREEN; // Vibrant green waveform color

    private final List<SoundTone> currentSoundTones;
    private final JList<String> soundToneList;
    private final DefaultListModel<String> soundToneListModel;
    private final JPanel soundToneParametersPanel;
    private final JPanel envelopeAndFilterPanel;

    private final JSlider durationSlider;
    private final JSlider offsetSlider;
    private final JSlider delayTimeSlider;
    private final JSlider delayDecaySlider;

    private final JSlider oscVolumeSlider;
    private final JSlider oscPitchSlider;
    private final JSlider oscDelaySlider;

    private final JComboBox<String> pitchEnvelopeFormComboBox;
    private final JTextField pitchEnvelopeStartField;
    private final JTextField pitchEnvelopeEndField;

    private final JComboBox<String> volumeEnvelopeFormComboBox;
    private final JTextField volumeEnvelopeStartField;
    private final JTextField volumeEnvelopeEndField;

    private final JCheckBox pitchModifierEnabledCheckbox;
    private final JComboBox<String> pitchModifierFormComboBox;
    private final JTextField pitchModifierStartField;
    private final JTextField pitchModifierEndField;

    private final JCheckBox pitchModifierAmplitudeEnabledCheckbox;
    private final JComboBox<String> pitchModifierAmplitudeFormComboBox;
    private final JTextField pitchModifierAmplitudeStartField;
    private final JTextField pitchModifierAmplitudeEndField;

    private final JCheckBox volumeMultiplierEnabledCheckbox;
    private final JComboBox<String> volumeMultiplierFormComboBox;
    private final JTextField volumeMultiplierStartField;
    private final JTextField volumeMultiplierEndField;

    private final JCheckBox volumeMultiplierAmplitudeEnabledCheckbox;
    private final JComboBox<String> volumeMultiplierAmplitudeFormComboBox;
    private final JTextField volumeMultiplierAmplitudeStartField;
    private final JTextField volumeMultiplierAmplitudeEndField;

    private final JCheckBox releaseEnabledCheckbox;
    private final JComboBox<String> releaseFormComboBox;
    private final JTextField releaseStartField;
    private final JTextField releaseEndField;

    private final JCheckBox attackEnabledCheckbox;
    private final JComboBox<String> attackFormComboBox;
    private final JTextField attackStartField;
    private final JTextField attackEndField;

    private final JCheckBox filterEnabledCheckbox;
    private final JSlider filterUnity0Slider;
    private final JSlider filterUnity1Slider;

    private final JSlider filterPairs0Slider; // For pairs[0]
    private final JSlider filterPairs1Slider; // For pairs[1]
    private final JPanel filterPairsPanel; // Panel to hold dynamic filter pair controls

    private Clip audioClip;
    private final ListSelectionListener listSelectionListener; // Store the ListSelectionListener

    private SoundTone currentlyEditedTone;

    private EnvelopeGraphPanel pitchEnvelopeGraphPanel;
    private EnvelopeGraphPanel volumeEnvelopeGraphPanel;
    private EnvelopeGraphPanel pitchModifierEnvelopeGraphPanel;
    private EnvelopeGraphPanel pitchModifierAmplitudeEnvelopeGraphPanel;
    private EnvelopeGraphPanel volumeMultiplierEnvelopeGraphPanel;
    private EnvelopeGraphPanel volumeMultiplierAmplitudeEnvelopeGraphPanel;
    private EnvelopeGraphPanel releaseEnvelopeGraphPanel;
    private EnvelopeGraphPanel attackEnvelopeGraphPanel;

    private final WaveformDisplayPanel waveformDisplayPanel;
    private FilterGraphPanel filterGraphPanel;

    private final Timer debounceTimer;

    private final JComboBox<String> outputBitDepthComboBox;
    private int outputBitDepth = 8;

    private final JList<String> oscillatorJList;
    private final int MAX_OSCILLATORS = 5;
    private int currentlySelectedOscillatorIndex = 0;


    public SoundEffectEditor() {
        setTitle("JagFX Editor - RuneScape Sound Effects v1");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setSize(1400, 900);
        setLocationRelativeTo(null);
        setLayout(new BorderLayout(5, 5));
        getContentPane().setBackground(BG_DEFAULT_FRAME);

        currentSoundTones = new ArrayList<>();
        currentSoundTones.add(new SoundTone());

        int DEBOUNCE_DELAY_MS = 100;
        debounceTimer = new Timer(DEBOUNCE_DELAY_MS, (ActionEvent _) -> {
            if (currentlyEditedTone != null) {
                updateWaveformDisplay();
                pitchEnvelopeGraphPanel.repaint();
                volumeEnvelopeGraphPanel.repaint();
                pitchModifierEnvelopeGraphPanel.repaint();
                pitchModifierAmplitudeEnvelopeGraphPanel.repaint();
                volumeMultiplierEnvelopeGraphPanel.repaint();
                volumeMultiplierAmplitudeEnvelopeGraphPanel.repaint();
                releaseEnvelopeGraphPanel.repaint();
                attackEnvelopeGraphPanel.repaint();
                filterGraphPanel.repaint();
            }
        });
        debounceTimer.setRepeats(false);

        JPanel topBar = new JPanel(new FlowLayout(FlowLayout.LEFT, 10, 10));
        topBar.setBackground(PANEL_BG);
        topBar.setBorder(new EmptyBorder(5, 10, 5, 10));

        topBar.add(createLabel("Sound Effect Controls:"));
        JButton playButton = createButton("Play Full Effect", BUTTON_GREEN);
        playButton.addActionListener(_ -> playSoundEffect());
        topBar.add(playButton);

        JButton playSelectedToneButton = createButton("Play Selected Tone", BUTTON_BLUE);
        playSelectedToneButton.addActionListener(_ -> playSelectedSoundTone());
        topBar.add(playSelectedToneButton);

        JButton stopButton = createButton("Stop Sound", BUTTON_RED);
        stopButton.addActionListener(_ -> stopSoundEffect());
        topBar.add(stopButton);

        JButton addToneButton = createButton("Add Tone", BUTTON_BLUE);
        addToneButton.addActionListener(_ -> addSoundTone());
        topBar.add(addToneButton);

        JButton removeToneButton = createButton("Remove Selected Tone", BUTTON_RED);
        removeToneButton.addActionListener(_ -> removeSelectedSoundTone());
        topBar.add(removeToneButton);

        JButton saveJagFXButton = createButton("Save as JagFX...", BUTTON_BLUE);
        saveJagFXButton.addActionListener(_ -> saveJagFXFile());
        topBar.add(saveJagFXButton);

        JButton importJagFXButton = createButton("Import JagFX...", BUTTON_BLUE);
        importJagFXButton.addActionListener(_ -> importJagFXFile());
        topBar.add(importJagFXButton);

        JButton exportWavButton = createButton("Export to WAV...", BUTTON_BLUE);
        exportWavButton.addActionListener(_ -> exportToWavFile());
        topBar.add(exportWavButton);

        topBar.add(createLabel("Output Bit Depth:"));
        outputBitDepthComboBox = createComboBox(new String[]{"8-bit", "16-bit"});
        outputBitDepthComboBox.setSelectedIndex(0); // Default to 8-bit
        outputBitDepthComboBox.addActionListener(_ -> {
            String selected = (String) outputBitDepthComboBox.getSelectedItem();
            if (selected != null) {
                outputBitDepth = selected.equals("8-bit") ? 8 : 16;
            }
            updateWaveformDisplay(); // Update waveform to reflect new bit depth scaling if needed
        });
        topBar.add(outputBitDepthComboBox);


        add(topBar, BorderLayout.NORTH);

        // --- Main Content Area (JSplitPane) ---
        JSplitPane mainSplitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
        mainSplitPane.setDividerSize(5);
        mainSplitPane.setResizeWeight(0.3); // More space for controls
        mainSplitPane.setBackground(BG_DEFAULT_FRAME);

        // Left Panel: Sound Tone List and Parameters
        JPanel leftControlPanel = new JPanel();
        leftControlPanel.setBackground(BG_DEFAULT_FRAME);
        leftControlPanel.setLayout(new BorderLayout(0, 5));
        leftControlPanel.setBorder(new EmptyBorder(5, 5, 5, 5));

        // Tone List
        JPanel toneListPanel = new JPanel(new BorderLayout(0, 5));
        toneListPanel.setBackground(PANEL_BG);
        toneListPanel.setBorder(BorderFactory.createCompoundBorder(
                new EmptyBorder(5, 5, 5, 5),
                new LineBorder(BORDER_COLOR, 1)
        ));
        toneListPanel.add(createTitleLabel("Sound Tones"), BorderLayout.NORTH);
        soundToneListModel = new DefaultListModel<>();
        soundToneList = new JList<>(soundToneListModel);
        soundToneList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        soundToneList.setBackground(PANEL_BG);
        soundToneList.setForeground(TEXT_COLOR_LIGHT);
        soundToneList.setFont(soundToneList.getFont().deriveFont(Font.PLAIN, 12f));
        soundToneList.setBorder(BorderFactory.createEmptyBorder());
        listSelectionListener = e -> {
            if (!e.getValueIsAdjusting()) {
                updateSoundToneEditorPanel();
            }
        };
        soundToneList.addListSelectionListener(listSelectionListener);
        JScrollPane toneListScrollPane = new JScrollPane(soundToneList);
        toneListScrollPane.setBorder(BorderFactory.createEmptyBorder());
        toneListPanel.add(toneListScrollPane, BorderLayout.CENTER);
        leftControlPanel.add(toneListPanel, BorderLayout.NORTH);


        JSplitPane centerRightSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
        centerRightSplitPane.setDividerSize(5);
        centerRightSplitPane.setResizeWeight(0.5);
        centerRightSplitPane.setBackground(BG_DEFAULT_FRAME);

        soundToneParametersPanel = new JPanel();
        soundToneParametersPanel.setBackground(PANEL_BG);
        soundToneParametersPanel.setLayout(new BorderLayout(0, 5));
        soundToneParametersPanel.setBorder(BorderFactory.createCompoundBorder(
                new EmptyBorder(5, 5, 5, 5),
                new LineBorder(BORDER_COLOR, 1)
        ));
        soundToneParametersPanel.add(createTitleLabel("Selected Sound Tone Parameters"), BorderLayout.NORTH);

        JPanel toneParamsContent = new JPanel(new GridBagLayout());
        toneParamsContent.setBackground(PANEL_BG);
        toneParamsContent.setBorder(new EmptyBorder(5, 5, 5, 5));
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(2, 5, 2, 5);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.WEST;

        gbc.gridx = 0; gbc.gridy = 0; toneParamsContent.add(createLabel("Duration (ms):"), gbc);
        gbc.gridx = 1; gbc.gridy = 0; gbc.weightx = 1.0; durationSlider = createSlider(0, 5000, 500); toneParamsContent.add(durationSlider, gbc);
        gbc.weightx = 0;

        gbc.gridx = 0; gbc.gridy = 1; toneParamsContent.add(createLabel("Offset (ms):"), gbc);
        gbc.gridx = 1; gbc.gridy = 1; gbc.weightx = 1.0; offsetSlider = createSlider(0, 2000, 0); toneParamsContent.add(offsetSlider, gbc);
        gbc.weightx = 0;

        gbc.gridx = 0; gbc.gridy = 2; toneParamsContent.add(createLabel("Delay Time (ms):"), gbc);
        gbc.gridx = 1; gbc.gridy = 2; gbc.weightx = 1.0; delayTimeSlider = createSlider(0, 1000, 0); toneParamsContent.add(delayTimeSlider, gbc);
        gbc.weightx = 0;

        gbc.gridx = 0; gbc.gridy = 3; toneParamsContent.add(createLabel("Delay Decay (%):"), gbc);
        gbc.gridx = 1; gbc.gridy = 3; gbc.weightx = 1.0; delayDecaySlider = createSlider(0, 100, 100); toneParamsContent.add(delayDecaySlider, gbc);
        gbc.weightx = 0;

        gbc.gridy = 4;
        gbc.gridwidth = 2;
        toneParamsContent.add(createTitleLabel("Oscillator Parameters"), gbc);
        gbc.gridwidth = 1;

        DefaultListModel<String> oscillatorListModel = new DefaultListModel<>();
        for (int i = 0; i < MAX_OSCILLATORS; i++) {
            oscillatorListModel.addElement("Oscillator " + (i + 1));
        }
        oscillatorJList = new JList<>(oscillatorListModel);
        oscillatorJList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        oscillatorJList.setBackground(PANEL_BG);
        oscillatorJList.setForeground(TEXT_COLOR_LIGHT);
        oscillatorJList.setFont(oscillatorJList.getFont().deriveFont(Font.PLAIN, 11f));
        oscillatorJList.setBorder(BorderFactory.createEmptyBorder());
        oscillatorJList.setFixedCellHeight(20);
        JScrollPane oscListScrollPane = new JScrollPane(oscillatorJList);
        oscListScrollPane.setPreferredSize(new Dimension(100, 120));
        oscListScrollPane.setBorder(BorderFactory.createLineBorder(BORDER_COLOR));

        oscillatorJList.addListSelectionListener(e -> {
            if (!e.getValueIsAdjusting()) {
                int selectedIndex = oscillatorJList.getSelectedIndex();
                if (selectedIndex != -1) {
                    currentlySelectedOscillatorIndex = selectedIndex;
                    updateOscillatorParametersPanel();
                }
            }
        });

        JPanel oscillatorListPanel = new JPanel(new BorderLayout());
        oscillatorListPanel.setBackground(PANEL_BG);
        oscillatorListPanel.add(oscListScrollPane, BorderLayout.CENTER);

        gbc.gridx = 0; gbc.gridy = 5; gbc.gridheight = 3; gbc.fill = GridBagConstraints.BOTH;
        toneParamsContent.add(oscillatorListPanel, gbc);
        gbc.gridheight = 1;

        JPanel currentOscillatorParametersPanel = new JPanel(new GridBagLayout());
        currentOscillatorParametersPanel.setBackground(PANEL_BG);
        currentOscillatorParametersPanel.setBorder(BorderFactory.createTitledBorder(new LineBorder(BORDER_COLOR), "Selected Oscillator"));
        GridBagConstraints oscGbc = new GridBagConstraints();
        oscGbc.insets = new Insets(2, 5, 2, 5);
        oscGbc.fill = GridBagConstraints.HORIZONTAL;
        oscGbc.anchor = GridBagConstraints.WEST;

        oscGbc.gridx = 0; oscGbc.gridy = 0; currentOscillatorParametersPanel.add(createLabel("Volume:"), oscGbc);
        oscGbc.gridx = 1; oscGbc.gridy = 0; oscGbc.weightx = 1.0; oscVolumeSlider = createSlider(0, 65535, 0); currentOscillatorParametersPanel.add(oscVolumeSlider, oscGbc); // Changed max to 65535
        oscGbc.weightx = 0;

        oscGbc.gridx = 0; oscGbc.gridy = 1; currentOscillatorParametersPanel.add(createLabel("Pitch:"), oscGbc);
        oscGbc.gridx = 1; oscGbc.gridy = 1; oscGbc.weightx = 1.0; oscPitchSlider = createSlider(-100, 100, 0); currentOscillatorParametersPanel.add(oscPitchSlider, oscGbc);
        oscGbc.weightx = 0;

        oscGbc.gridx = 0; oscGbc.gridy = 2; currentOscillatorParametersPanel.add(createLabel("Delay:"), oscGbc);
        oscGbc.gridx = 1; oscGbc.gridy = 2; oscGbc.weightx = 1.0; oscDelaySlider = createSlider(0, 9, 0); currentOscillatorParametersPanel.add(oscDelaySlider, oscGbc); // Delay type 0-9
        oscGbc.weightx = 0;

        gbc.gridx = 1; gbc.gridy = 5; gbc.gridheight = 3; gbc.weightx = 1.0; gbc.fill = GridBagConstraints.BOTH;
        toneParamsContent.add(currentOscillatorParametersPanel, gbc);
        gbc.gridheight = 1; gbc.weightx = 0; // Reset gridheight and weightx

        soundToneParametersPanel.add(new JScrollPane(toneParamsContent), BorderLayout.CENTER);
        centerRightSplitPane.setTopComponent(soundToneParametersPanel);

        // Right Panel: Envelope and Filter Parameters
        envelopeAndFilterPanel = new JPanel();
        envelopeAndFilterPanel.setBackground(PANEL_BG);
        envelopeAndFilterPanel.setLayout(new BoxLayout(envelopeAndFilterPanel, BoxLayout.Y_AXIS));
        envelopeAndFilterPanel.setBorder(BorderFactory.createCompoundBorder(
                new EmptyBorder(5, 5, 5, 5),
                new LineBorder(BORDER_COLOR, 1)
        ));

        // Pitch Envelope
        envelopeAndFilterPanel.add(createEnvelopePanel("Pitch Envelope",
                pitchEnvelopeFormComboBox = createComboBox(new String[]{"Square", "Sine", "Triangle", "Noise", "Off"}), // Added "Off"
                pitchEnvelopeStartField = createTextField("0"),
                pitchEnvelopeEndField = createTextField("0"),
                null // No checkbox for main envelopes
        ));

        // Volume Envelope
        envelopeAndFilterPanel.add(createEnvelopePanel("Volume Envelope",
                volumeEnvelopeFormComboBox = createComboBox(new String[]{"Square", "Sine", "Triangle", "Noise", "Off"}), // Added "Off"
                volumeEnvelopeStartField = createTextField("0"),
                volumeEnvelopeEndField = createTextField("0"),
                null // No checkbox for main envelopes
        ));

        // NEW: Pitch Modifier Envelope
        envelopeAndFilterPanel.add(createEnvelopePanel("Pitch Modifier Envelope",
                pitchModifierFormComboBox = createComboBox(new String[]{"Square", "Sine", "Triangle", "Noise", "Off"}),
                pitchModifierStartField = createTextField("0"),
                pitchModifierEndField = createTextField("0"),
                pitchModifierEnabledCheckbox = new JCheckBox("Enable Pitch Modifier")
        ));

        // NEW: Pitch Modifier Amplitude Envelope
        envelopeAndFilterPanel.add(createEnvelopePanel("Pitch Modifier Amplitude Envelope",
                pitchModifierAmplitudeFormComboBox = createComboBox(new String[]{"Square", "Sine", "Triangle", "Noise", "Off"}),
                pitchModifierAmplitudeStartField = createTextField("0"),
                pitchModifierAmplitudeEndField = createTextField("0"),
                pitchModifierAmplitudeEnabledCheckbox = new JCheckBox("Enable Pitch Mod Amplitude")
        ));

        // NEW: Volume Multiplier Envelope
        envelopeAndFilterPanel.add(createEnvelopePanel("Volume Multiplier Envelope",
                volumeMultiplierFormComboBox = createComboBox(new String[]{"Square", "Sine", "Triangle", "Noise", "Off"}),
                volumeMultiplierStartField = createTextField("0"),
                volumeMultiplierEndField = createTextField("0"),
                volumeMultiplierEnabledCheckbox = new JCheckBox("Enable Volume Multiplier")
        ));

        // NEW: Volume Multiplier Amplitude Envelope
        envelopeAndFilterPanel.add(createEnvelopePanel("Volume Multiplier Amplitude Envelope",
                volumeMultiplierAmplitudeFormComboBox = createComboBox(new String[]{"Square", "Sine", "Triangle", "Noise", "Off"}),
                volumeMultiplierAmplitudeStartField = createTextField("0"),
                volumeMultiplierAmplitudeEndField = createTextField("0"),
                volumeMultiplierAmplitudeEnabledCheckbox = new JCheckBox("Enable Volume Mul Amplitude")
        ));

        // NEW: Release Envelope
        envelopeAndFilterPanel.add(createEnvelopePanel("Release Envelope",
                releaseFormComboBox = createComboBox(new String[]{"Square", "Sine", "Triangle", "Noise", "Off"}),
                releaseStartField = createTextField("0"),
                releaseEndField = createTextField("0"),
                releaseEnabledCheckbox = new JCheckBox("Enable Release")
        ));

        // NEW: Attack Envelope
        envelopeAndFilterPanel.add(createEnvelopePanel("Attack Envelope",
                attackFormComboBox = createComboBox(new String[]{"Square", "Sine", "Triangle", "Noise", "Off"}),
                attackStartField = createTextField("0"),
                attackEndField = createTextField("0"),
                attackEnabledCheckbox = new JCheckBox("Enable Attack")
        ));


        // Sound Filter
        JPanel filterPanel = new JPanel(new GridBagLayout());
        filterPanel.setBackground(PANEL_BG);
        filterPanel.setBorder(BorderFactory.createTitledBorder(new LineBorder(BORDER_COLOR), "Sound Filter"));
        int filterRow = 0;

        gbc.gridx = 0; gbc.gridy = filterRow; filterEnabledCheckbox = new JCheckBox("Enable Filter");
        filterEnabledCheckbox.setBackground(PANEL_BG);
        filterEnabledCheckbox.setForeground(TEXT_COLOR_LIGHT);
        filterEnabledCheckbox.setFont(filterEnabledCheckbox.getFont().deriveFont(Font.PLAIN, 11f));
        gbc.gridwidth = 2; filterPanel.add(filterEnabledCheckbox, gbc);
        gbc.gridwidth = 1; filterRow++;

        gbc.gridx = 0; gbc.gridy = filterRow; filterPanel.add(createLabel("Unity 0:"), gbc);
        gbc.gridx = 1; gbc.gridy = filterRow++; filterUnity0Slider = createSlider(0, 65535, 0); filterPanel.add(filterUnity0Slider, gbc);

        gbc.gridx = 0; gbc.gridy = filterRow; filterPanel.add(createLabel("Unity 1:"), gbc);
        gbc.gridx = 1; gbc.gridy = filterRow++; filterUnity1Slider = createSlider(0, 65535, 0); filterPanel.add(filterUnity1Slider, gbc);

        // NEW: Filter Pairs
        gbc.gridx = 0; gbc.gridy = filterRow; filterPanel.add(createLabel("Pairs 0:"), gbc);
        gbc.gridx = 1; gbc.gridy = filterRow++; filterPairs0Slider = createSlider(0, 4, 0); filterPanel.add(filterPairs0Slider, gbc); // Max 4 pairs
        filterPairs0Slider.setMajorTickSpacing(1); // Set major tick spacing to 1
        filterPairs0Slider.setMinorTickSpacing(1); // Set minor tick spacing to 1
        filterPairs0Slider.setPaintLabels(true); // Show labels for ticks

        gbc.gridx = 0; gbc.gridy = filterRow; filterPanel.add(createLabel("Pairs 1:"), gbc);
        gbc.gridx = 1; gbc.gridy = filterRow++; filterPairs1Slider = createSlider(0, 4, 0); filterPanel.add(filterPairs1Slider, gbc); // Max 4 pairs
        filterPairs1Slider.setMajorTickSpacing(1); // Set major tick spacing to 1
        filterPairs1Slider.setMinorTickSpacing(1); // Set minor tick spacing to 1
        filterPairs1Slider.setPaintLabels(true); // Show labels for ticks

        filterPairsPanel = new JPanel(new GridBagLayout());
        filterPairsPanel.setBackground(PANEL_BG);
        filterPairsPanel.setBorder(BorderFactory.createTitledBorder(new LineBorder(BORDER_COLOR), "Filter Coefficients"));
        gbc.gridx = 0; gbc.gridy = filterRow;
        gbc.gridwidth = 2; gbc.weighty = 1.0; gbc.fill = GridBagConstraints.BOTH;
        filterPanel.add(new JScrollPane(filterPairsPanel), gbc); // Wrap in scroll pane
        gbc.gridwidth = 1; gbc.weighty = 0; // Reset

        envelopeAndFilterPanel.add(filterPanel);

        centerRightSplitPane.setBottomComponent(new JScrollPane(envelopeAndFilterPanel)); // Wrap in scroll pane
        leftControlPanel.add(centerRightSplitPane, BorderLayout.CENTER); // Add the split pane to the left control panel
        mainSplitPane.setLeftComponent(leftControlPanel);

        // Right Panel: Graphs and Waveform
        JPanel rightGraphPanel = new JPanel();
        rightGraphPanel.setBackground(BG_DEFAULT_FRAME);
        rightGraphPanel.setLayout(new GridLayout(4, 2, 5, 5)); // 4 rows, 2 columns for graphs (more space)
        rightGraphPanel.setBorder(new EmptyBorder(5, 5, 5, 5));

        // Initialize graph panels
        pitchEnvelopeGraphPanel = new EnvelopeGraphPanel("Pitch Envelope", GRAPH_BG, GRAPH_GRID, GRAPH_LINE);
        volumeEnvelopeGraphPanel = new EnvelopeGraphPanel("Volume Envelope", GRAPH_BG, GRAPH_GRID, GRAPH_LINE);
        pitchModifierEnvelopeGraphPanel = new EnvelopeGraphPanel("Pitch Modulator", GRAPH_BG, GRAPH_GRID, GRAPH_LINE);
        pitchModifierAmplitudeEnvelopeGraphPanel = new EnvelopeGraphPanel("Pitch Mod Amplitude", GRAPH_BG, GRAPH_GRID, GRAPH_LINE);
        volumeMultiplierEnvelopeGraphPanel = new EnvelopeGraphPanel("Volume Multiplier", GRAPH_BG, GRAPH_GRID, GRAPH_LINE);
        volumeMultiplierAmplitudeEnvelopeGraphPanel = new EnvelopeGraphPanel("Volume Mul Amplitude", GRAPH_BG, GRAPH_GRID, GRAPH_LINE);
        releaseEnvelopeGraphPanel = new EnvelopeGraphPanel("Release Envelope", GRAPH_BG, GRAPH_GRID, GRAPH_LINE);
        attackEnvelopeGraphPanel = new EnvelopeGraphPanel("Attack Envelope", GRAPH_BG, GRAPH_GRID, GRAPH_LINE);

        waveformDisplayPanel = new WaveformDisplayPanel("Waveform Output", WAVEFORM_BG, WAVEFORM_COLOR);
        filterGraphPanel = new FilterGraphPanel("Filter Frequency Response", GRAPH_BG, GRAPH_GRID, GRAPH_LINE);

        rightGraphPanel.add(pitchEnvelopeGraphPanel);
        rightGraphPanel.add(volumeEnvelopeGraphPanel);
        rightGraphPanel.add(pitchModifierEnvelopeGraphPanel);
        rightGraphPanel.add(pitchModifierAmplitudeEnvelopeGraphPanel);
        rightGraphPanel.add(volumeMultiplierEnvelopeGraphPanel);
        rightGraphPanel.add(volumeMultiplierAmplitudeEnvelopeGraphPanel);
        rightGraphPanel.add(releaseEnvelopeGraphPanel);
        rightGraphPanel.add(attackEnvelopeGraphPanel);
        rightGraphPanel.add(filterGraphPanel);
        rightGraphPanel.add(waveformDisplayPanel);


        mainSplitPane.setRightComponent(rightGraphPanel);
        add(mainSplitPane, BorderLayout.CENTER);

        // Initialize the tone list and select the first tone
        updateToneList();
        if (!currentSoundTones.isEmpty()) {
            soundToneList.setSelectedIndex(0);
        }

        // Initialize oscillator selection
        oscillatorJList.setSelectedIndex(0);

        // Add listeners once during initialization
        addGlobalParameterListeners();

        setVisible(true);
    }

    private JLabel createLabel(String text) {
        JLabel label = new JLabel(text);
        label.setForeground(TEXT_COLOR_LIGHT);
        label.setFont(label.getFont().deriveFont(Font.PLAIN, 12f));
        return label;
    }

    private JLabel createTitleLabel(String text) {
        JLabel label = new JLabel(text);
        label.setForeground(HIGHLIGHT_COLOR_GREEN);
        label.setFont(label.getFont().deriveFont(Font.BOLD, 14f));
        label.setBorder(new EmptyBorder(5, 5, 5, 5));
        return label;
    }

    private JButton createButton(String text, Color bgColor) {
        JButton button = new JButton(text);
        button.setBackground(bgColor);
        button.setForeground(Color.WHITE);
        button.setFont(button.getFont().deriveFont(Font.BOLD, 11f));
        button.setFocusPainted(false);
        button.setBorder(BorderFactory.createCompoundBorder(
                new LineBorder(BORDER_COLOR, 1),
                new EmptyBorder(5, 10, 5, 10)
        ));
        button.setOpaque(true);
        button.setCursor(new Cursor(Cursor.HAND_CURSOR));

        if (bgColor.equals(BUTTON_BG)) {
            button.setForeground(TEXT_COLOR_LIGHT);
        }
        return button;
    }

    private JComboBox<String> createComboBox(String[] items) {
        JComboBox<String> comboBox = new JComboBox<>(items);
        comboBox.setBackground(BUTTON_BG);
        comboBox.setForeground(TEXT_COLOR_LIGHT);
        comboBox.setFont(comboBox.getFont().deriveFont(Font.PLAIN, 12f));
        comboBox.setBorder(BorderFactory.createLineBorder(BORDER_COLOR));
        comboBox.setRenderer(new DefaultListCellRenderer() {
            @Override
            public Component getListCellRendererComponent(JList<?> list, Object value, int index, boolean isSelected, boolean cellHasFocus) {
                JLabel label = (JLabel) super.getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);
                label.setBackground(isSelected ? UIManager.getColor("ComboBox.selectionBackground") : BUTTON_BG);
                label.setForeground(isSelected ? UIManager.getColor("ComboBox.selectionForeground") : TEXT_COLOR_LIGHT);
                label.setBorder(new EmptyBorder(2, 5, 2, 5));
                return label;
            }
        });
        return comboBox;
    }

    private JTextField createTextField(String text) {
        JTextField textField = new JTextField(text);
        textField.setBackground(BUTTON_BG);
        textField.setForeground(TEXT_COLOR_LIGHT);
        textField.setFont(textField.getFont().deriveFont(Font.PLAIN, 12f));
        textField.setBorder(BorderFactory.createLineBorder(BORDER_COLOR));
        textField.setPreferredSize(new Dimension(80, 25));
        return textField;
    }

    private JSlider createSlider(int min, int max, int value) {
        JSlider slider = new JSlider(min, max, value);
        slider.setBackground(PANEL_BG);
        slider.setForeground(TEXT_COLOR_LIGHT);
        slider.setPaintTicks(true);
        slider.setPaintLabels(true);
        slider.setMajorTickSpacing((max - min) / 4);
        slider.setMinorTickSpacing((max - min) / 20);
        slider.setFocusable(false);

        slider.setUI(new BasicSliderUI(slider) {
            @Override
            public void paintTrack(Graphics g) {
                Graphics2D g2d = (Graphics2D) g;
                Rectangle trackBounds = trackRect;
                g2d.setColor(SLIDER_TRACK_COLOR);
                g2d.fillRoundRect(trackBounds.x, trackBounds.y + trackBounds.height / 2 - 2, trackBounds.width, 4, 4, 4);
            }

            @Override
            public void paintThumb(Graphics g) {
                Graphics2D g2d = (Graphics2D) g;
                Rectangle thumbBounds = thumbRect;
                g2d.setColor(SLIDER_THUMB_COLOR);
                g2d.fillOval(thumbBounds.x, thumbBounds.y, thumbBounds.width, thumbBounds.height);
                g2d.setColor(BORDER_COLOR);
                g2d.drawOval(thumbBounds.x, thumbBounds.y, thumbBounds.width, thumbBounds.height);
            }
        });
        return slider;
    }

    // Helper to create an envelope panel
    private JPanel createEnvelopePanel(String title, JComboBox<String> formComboBox, JTextField startField, JTextField endField, JCheckBox enabledCheckbox) {
        JPanel panel = new JPanel(new GridBagLayout());
        panel.setBackground(PANEL_BG);
        panel.setBorder(BorderFactory.createTitledBorder(
                new LineBorder(BORDER_COLOR),
                title,
                TitledBorder.LEFT,
                TitledBorder.TOP,
                new Font("Arial", Font.BOLD, 12),
                HIGHLIGHT_COLOR_GREEN // Changed to green highlight
        ));
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(2, 5, 2, 5);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.WEST;
        int row = 0;

        if (enabledCheckbox != null) {
            gbc.gridx = 0; gbc.gridy = row; gbc.gridwidth = 2;
            enabledCheckbox.setBackground(PANEL_BG);
            enabledCheckbox.setForeground(TEXT_COLOR_LIGHT); // Changed to light text
            enabledCheckbox.setFont(enabledCheckbox.getFont().deriveFont(Font.PLAIN, 11f));
            panel.add(enabledCheckbox, gbc);
            gbc.gridwidth = 1; row++;
        }

        gbc.gridx = 0; gbc.gridy = row; panel.add(createLabel("Form:"), gbc);
        gbc.gridx = 1; gbc.gridy = row++; panel.add(formComboBox, gbc);

        gbc.gridx = 0; gbc.gridy = row; panel.add(createLabel("Start:"), gbc);
        gbc.gridx = 1; gbc.gridy = row++; panel.add(startField, gbc);

        gbc.gridx = 0; gbc.gridy = row; panel.add(createLabel("End:"), gbc);
        gbc.gridx = 1; gbc.gridy = row;
        panel.add(endField, gbc);

        return panel;
    }

    private void addSoundTone() {
        if (currentSoundTones.size() < 10) {
            currentSoundTones.add(new SoundTone());
            updateToneList(); // This now handles setting selection and calling updateSoundToneEditorPanel()
        } else {
            JOptionPane.showMessageDialog(this, "Maximum of 10 Sound Tones reached.", "Limit Reached", JOptionPane.WARNING_MESSAGE);
        }
    }

    private void removeSelectedSoundTone() {
        int selectedIndex = soundToneList.getSelectedIndex();
        if (selectedIndex != -1 && currentSoundTones.size() > 1) { // Ensure at least one tone remains
            currentSoundTones.remove(selectedIndex);
            updateToneList(); // This now handles setting selection and calling updateSoundToneEditorPanel()
        } else if (currentSoundTones.size() == 1) {
            JOptionPane.showMessageDialog(this, "Cannot remove the last Sound Tone.", "Minimum Tones", JOptionPane.WARNING_MESSAGE);
        } else {
            JOptionPane.showMessageDialog(this, "No Sound Tone selected to remove.", "No Selection", JOptionPane.WARNING_MESSAGE);
        }
    }

    private void updateToneList() {
        new Thread(() -> {
            // Remove the listener to prevent it from firing during model updates
            soundToneList.removeListSelectionListener(listSelectionListener);

            soundToneListModel.clear();
            for (int i = 0; i < currentSoundTones.size(); i++) {
                soundToneListModel.addElement("Tone " + (i + 1));
            }

            // Re-add the listener
            soundToneList.addListSelectionListener(listSelectionListener);

            // Set the selection. This will trigger the listSelectionListener.
            if (!currentSoundTones.isEmpty()) {
                int newSelectionIndex = soundToneList.getSelectedIndex();
                if (newSelectionIndex == -1 || newSelectionIndex >= currentSoundTones.size()) {
                    newSelectionIndex = 0; // Default to selecting the first tone
                }
                soundToneList.setSelectedIndex(newSelectionIndex);
            } else {
                // If no tones left, ensure the panel is cleared.
                // The listener won't fire if the list becomes empty.
                clearSoundToneEditorPanel();
            }
        }).start();
    }

    private void updateSoundToneEditorPanel() {
        int selectedIndex = soundToneList.getSelectedIndex();
        if (selectedIndex != -1) {
            currentlyEditedTone = currentSoundTones.get(selectedIndex); // Set the currently edited tone

            // Update main SoundTone parameters
            durationSlider.setValue(currentlyEditedTone.duration);
            offsetSlider.setValue(currentlyEditedTone.offset);
            delayTimeSlider.setValue(currentlyEditedTone.delayTime);
            delayDecaySlider.setValue(currentlyEditedTone.delayDecay);

            updateOscillatorParametersPanel();

            updateEnvelopeUI(currentlyEditedTone.pitch, pitchEnvelopeFormComboBox, pitchEnvelopeStartField, pitchEnvelopeEndField, pitchEnvelopeGraphPanel, null);
            updateEnvelopeUI(currentlyEditedTone.volume, volumeEnvelopeFormComboBox, volumeEnvelopeStartField, volumeEnvelopeEndField, volumeEnvelopeGraphPanel, null);
            updateEnvelopeUI(currentlyEditedTone.pitchModifier, pitchModifierFormComboBox, pitchModifierStartField, pitchModifierEndField, pitchModifierEnvelopeGraphPanel, pitchModifierEnabledCheckbox);
            updateEnvelopeUI(currentlyEditedTone.pitchModifierAmplitude, pitchModifierAmplitudeFormComboBox, pitchModifierAmplitudeStartField, pitchModifierAmplitudeEndField, pitchModifierAmplitudeEnvelopeGraphPanel, pitchModifierAmplitudeEnabledCheckbox);
            updateEnvelopeUI(currentlyEditedTone.volumeMultiplier, volumeMultiplierFormComboBox, volumeMultiplierStartField, volumeMultiplierEndField, volumeMultiplierEnvelopeGraphPanel, volumeMultiplierEnabledCheckbox);
            updateEnvelopeUI(currentlyEditedTone.volumeMultiplierAmplitude, volumeMultiplierAmplitudeFormComboBox, volumeMultiplierAmplitudeStartField, volumeMultiplierAmplitudeEndField, volumeMultiplierAmplitudeEnvelopeGraphPanel, volumeMultiplierAmplitudeEnabledCheckbox);
            updateEnvelopeUI(currentlyEditedTone.release, releaseFormComboBox, releaseStartField, releaseEndField, releaseEnvelopeGraphPanel, releaseEnabledCheckbox);
            updateEnvelopeUI(currentlyEditedTone.attack, attackFormComboBox, attackStartField, attackEndField, attackEnvelopeGraphPanel, attackEnabledCheckbox);

            if (currentlyEditedTone.filter != null) {
                boolean filterEnabled = currentlyEditedTone.filter.pairs[0] > 0 || currentlyEditedTone.filter.pairs[1] > 0;
                filterEnabledCheckbox.setSelected(filterEnabled);

                filterUnity0Slider.setValue(currentlyEditedTone.filter.unity[0]);
                filterUnity1Slider.setValue(currentlyEditedTone.filter.unity[1]);
                // Clamp the slider values to 4, even if the internal 'pairs' might be higher from a malformed file
                filterPairs0Slider.setValue(Math.min(currentlyEditedTone.filter.pairs[0], 4));
                filterPairs1Slider.setValue(Math.min(currentlyEditedTone.filter.pairs[1], 4));

                // Enable/disable based on filterEnabledCheckbox
                filterUnity0Slider.setEnabled(filterEnabled);
                filterUnity1Slider.setEnabled(filterEnabled);
                filterPairs0Slider.setEnabled(filterEnabled);
                filterPairs1Slider.setEnabled(filterEnabled);
                updateFilterPairsPanel(filterEnabled); // Update dynamic filter controls

                filterGraphPanel.setFilter(currentlyEditedTone.filter);
            } else {
                filterEnabledCheckbox.setSelected(false);
                filterUnity0Slider.setValue(0);
                filterUnity1Slider.setValue(0);
                filterPairs0Slider.setValue(0);
                filterPairs1Slider.setValue(0);
                filterUnity0Slider.setEnabled(false);
                filterUnity1Slider.setEnabled(false);
                filterPairs0Slider.setEnabled(false);
                filterPairs1Slider.setEnabled(false);
                updateFilterPairsPanel(false); // Clear dynamic filter controls
                filterGraphPanel.setFilter(null);
            }

        } else {
            clearSoundToneEditorPanel();
            currentlyEditedTone = null; // Clear the reference
        }
        // Revalidate and repaint to ensure UI updates are rendered
        soundToneParametersPanel.revalidate();
        soundToneParametersPanel.repaint();
        envelopeAndFilterPanel.revalidate();
        envelopeAndFilterPanel.repaint();
    }

    private void updateEnvelopeUI(SoundEnvelope envelope, JComboBox<String> formComboBox, JTextField startField, JTextField endField, EnvelopeGraphPanel graphPanel, JCheckBox enabledCheckbox) {
        if (enabledCheckbox != null) {
            if (envelope != null) {
                enabledCheckbox.setSelected(envelope.form != 0);
                int formIndex = envelope.form == 0 ? 4 : envelope.form - 1;
                if (formIndex >= 0 && formIndex < formComboBox.getItemCount()) {
                    formComboBox.setSelectedIndex(formIndex);
                } else {
                    formComboBox.setSelectedIndex(4); // Default to "Off" if invalid
                }
                startField.setText(String.valueOf(envelope.start));
                endField.setText(String.valueOf(envelope.end));
                graphPanel.setEnvelope(envelope);
                // Enable controls if checkbox is selected
                boolean isActive = enabledCheckbox.isSelected();
                formComboBox.setEnabled(isActive);
                startField.setEnabled(isActive);
                endField.setEnabled(isActive);
            } else { // If envelope object is null (not present in file)
                enabledCheckbox.setSelected(false);
                formComboBox.setSelectedIndex(4); // "Off"
                startField.setText("0");
                endField.setText("0");
                graphPanel.setEnvelope(null);
                // Disable all controls
                formComboBox.setEnabled(false);
                startField.setEnabled(false);
                endField.setEnabled(false);
            }
        } else { // For always-present envelopes (pitch, volume)
            if (envelope != null) {
                int formIndex = envelope.form == 0 ? 4 : envelope.form - 1;
                if (formIndex >= 0 && formIndex < formComboBox.getItemCount()) {
                    formComboBox.setSelectedIndex(formIndex);
                } else {
                    formComboBox.setSelectedIndex(4); // Default to "Off" if invalid
                }
                startField.setText(String.valueOf(envelope.start));
                endField.setText(String.valueOf(envelope.end));
                graphPanel.setEnvelope(envelope);
                // Always enabled for main envelopes
                formComboBox.setEnabled(true);
                startField.setEnabled(true);
                endField.setEnabled(true);
            } else {
                // This case should ideally not happen for pitch/volume envelopes
                // but as a fallback, clear and disable.
                formComboBox.setSelectedIndex(4); // "Off"
                startField.setText("0");
                endField.setText("0");
                graphPanel.setEnvelope(null);
                formComboBox.setEnabled(false);
                startField.setEnabled(false);
                endField.setEnabled(false);
            }
        }
    }


    private void clearSoundToneEditorPanel() {
        durationSlider.setValue(0);
        offsetSlider.setValue(0);
        delayTimeSlider.setValue(0);
        delayDecaySlider.setValue(0);

        // Clear oscillator parameters
        currentlySelectedOscillatorIndex = 0;
        updateOscillatorParametersPanel(); // This will set current osc params to 0

        // For always-present envelopes, pass a new default instance
        updateEnvelopeUI(new SoundEnvelope(), pitchEnvelopeFormComboBox, pitchEnvelopeStartField, pitchEnvelopeEndField, pitchEnvelopeGraphPanel, null);
        updateEnvelopeUI(new SoundEnvelope(), volumeEnvelopeFormComboBox, volumeEnvelopeStartField, volumeEnvelopeEndField, volumeEnvelopeGraphPanel, null);

        // For optional envelopes, pass null to clear and disable
        updateEnvelopeUI(null, pitchModifierFormComboBox, pitchModifierStartField, pitchModifierEndField, pitchModifierEnvelopeGraphPanel, pitchModifierEnabledCheckbox);
        updateEnvelopeUI(null, pitchModifierAmplitudeFormComboBox, pitchModifierAmplitudeStartField, pitchModifierAmplitudeEndField, pitchModifierAmplitudeEnvelopeGraphPanel, pitchModifierAmplitudeEnabledCheckbox);
        updateEnvelopeUI(null, volumeMultiplierFormComboBox, volumeMultiplierStartField, volumeMultiplierEndField, volumeMultiplierEnvelopeGraphPanel, volumeMultiplierEnabledCheckbox);
        updateEnvelopeUI(null, volumeMultiplierAmplitudeFormComboBox, volumeMultiplierAmplitudeStartField, volumeMultiplierAmplitudeEndField, volumeMultiplierAmplitudeEnvelopeGraphPanel, volumeMultiplierAmplitudeEnabledCheckbox);
        updateEnvelopeUI(null, releaseFormComboBox, releaseStartField, releaseEndField, releaseEnvelopeGraphPanel, releaseEnabledCheckbox);
        updateEnvelopeUI(null, attackFormComboBox, attackStartField, attackEndField, attackEnvelopeGraphPanel, attackEnabledCheckbox);

        filterEnabledCheckbox.setSelected(false);
        filterUnity0Slider.setValue(0);
        filterUnity1Slider.setValue(0);
        filterPairs0Slider.setValue(0);
        filterPairs1Slider.setValue(0);
        filterUnity0Slider.setEnabled(false);
        filterUnity1Slider.setEnabled(false);
        filterPairs0Slider.setEnabled(false);
        filterPairs1Slider.setEnabled(false);
        updateFilterPairsPanel(false); // Clear dynamic filter controls
        filterGraphPanel.setFilter(null);

        waveformDisplayPanel.setAudioData(new byte[0]); // Clear waveform data
    }

    // Update the UI for the currently selected oscillator
    private void updateOscillatorParametersPanel() {
        if (currentlyEditedTone != null && currentlySelectedOscillatorIndex >= 0 && currentlySelectedOscillatorIndex < MAX_OSCILLATORS) {
            oscVolumeSlider.setValue(currentlyEditedTone.oscillatorVolume[currentlySelectedOscillatorIndex]);
            oscPitchSlider.setValue(currentlyEditedTone.oscillatorPitch[currentlySelectedOscillatorIndex]);
            oscDelaySlider.setValue(currentlyEditedTone.oscillatorDelays[currentlySelectedOscillatorIndex]);
        } else {
            oscVolumeSlider.setValue(0);
            oscPitchSlider.setValue(0);
            oscDelaySlider.setValue(0);
        }
    }

    // Dynamically update the filter pairs panel
    private void updateFilterPairsPanel(boolean enabled) {
        filterPairsPanel.removeAll(); // Clear existing controls
        if (currentlyEditedTone == null || !enabled || currentlyEditedTone.filter == null) {
            filterPairsPanel.revalidate();
            filterPairsPanel.repaint();
            return;
        }

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(2, 5, 2, 5);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.WEST;

        // Direction 0 (Input)
        gbc.gridx = 0; gbc.gridy = 0; gbc.gridwidth = 4;
        filterPairsPanel.add(createLabel("Direction 0 (Input)"), gbc);
        gbc.gridwidth = 1;

        // Clamp the iteration count for UI elements to the array size (4)
        int numPairs0 = Math.min(currentlyEditedTone.filter.pairs[0], 4);
        for (int i = 0; i < numPairs0; i++) {
            final int pairIndex = i;
            gbc.gridy++;
            filterPairsPanel.add(createLabel("Pair " + (i + 1) + " (0):"), gbc);

            gbc.gridx = 1; filterPairsPanel.add(createLabel("Phase 0:"), gbc);
            JTextField phase0Field = createTextField(String.valueOf(currentlyEditedTone.filter.phases[0][0][i]));
            phase0Field.getDocument().addDocumentListener(new DocumentListener() {
                public void changedUpdate(DocumentEvent e) { updateFilterPhase(0, 0, pairIndex, phase0Field); }
                public void removeUpdate(DocumentEvent e) { updateFilterPhase(0, 0, pairIndex, phase0Field); }
                public void insertUpdate(DocumentEvent e) { updateFilterPhase(0, 0, pairIndex, phase0Field); }
            });
            gbc.gridx = 2; filterPairsPanel.add(phase0Field, gbc);

            gbc.gridx = 3; filterPairsPanel.add(createLabel("Mag 0:"), gbc);
            JTextField mag0Field = createTextField(String.valueOf(currentlyEditedTone.filter.magnitudes[0][0][i]));
            mag0Field.getDocument().addDocumentListener(new DocumentListener() {
                public void changedUpdate(DocumentEvent e) { updateFilterMagnitude(0, 0, pairIndex, mag0Field); }
                public void removeUpdate(DocumentEvent e) { updateFilterMagnitude(0, 0, pairIndex, mag0Field); }
                public void insertUpdate(DocumentEvent e) { updateFilterMagnitude(0, 0, pairIndex, mag0Field); }
            });
            gbc.gridx = 4; filterPairsPanel.add(mag0Field, gbc);

            gbc.gridx = 1; gbc.gridy++; filterPairsPanel.add(createLabel("Phase 1:"), gbc);
            JTextField phase1Field = createTextField(String.valueOf(currentlyEditedTone.filter.phases[0][1][i]));
            phase1Field.getDocument().addDocumentListener(new DocumentListener() {
                public void changedUpdate(DocumentEvent e) { updateFilterPhase(0, 1, pairIndex, phase1Field); }
                public void removeUpdate(DocumentEvent e) { updateFilterPhase(0, 1, pairIndex, phase1Field); }
                public void insertUpdate(DocumentEvent e) { updateFilterPhase(0, 1, pairIndex, phase1Field); }
            });
            gbc.gridx = 2; filterPairsPanel.add(phase1Field, gbc);

            gbc.gridx = 3; filterPairsPanel.add(createLabel("Mag 1:"), gbc);
            JTextField mag1Field = createTextField(String.valueOf(currentlyEditedTone.filter.magnitudes[0][1][i]));
            mag1Field.getDocument().addDocumentListener(new DocumentListener() {
                public void changedUpdate(DocumentEvent e) { updateFilterMagnitude(0, 1, pairIndex, mag1Field); }
                public void removeUpdate(DocumentEvent e) { updateFilterMagnitude(0, 1, pairIndex, mag1Field); }
                public void insertUpdate(DocumentEvent e) { updateFilterMagnitude(0, 1, pairIndex, mag1Field); }
            });
            gbc.gridx = 4; filterPairsPanel.add(mag1Field, gbc);
        }

        // Direction 1 (Output)
        gbc.gridx = 0; gbc.gridy++; gbc.gridwidth = 4;
        filterPairsPanel.add(createLabel("Direction 1 (Output)"), gbc);
        gbc.gridwidth = 1;

        // Clamp the iteration count for UI elements to the array size (4)
        int numPairs1 = Math.min(currentlyEditedTone.filter.pairs[1], 4);
        for (int i = 0; i < numPairs1; i++) {
            final int pairIndex = i;
            gbc.gridy++;
            filterPairsPanel.add(createLabel("Pair " + (i + 1) + " (1):"), gbc);

            gbc.gridx = 1; filterPairsPanel.add(createLabel("Phase 0:"), gbc);
            JTextField phase0Field = createTextField(String.valueOf(currentlyEditedTone.filter.phases[1][0][i]));
            phase0Field.getDocument().addDocumentListener(new DocumentListener() {
                public void changedUpdate(DocumentEvent e) { updateFilterPhase(1, 0, pairIndex, phase0Field); }
                public void removeUpdate(DocumentEvent e) { updateFilterPhase(1, 0, pairIndex, phase0Field); }
                public void insertUpdate(DocumentEvent e) { updateFilterPhase(1, 0, pairIndex, phase0Field); }
            });
            gbc.gridx = 2; filterPairsPanel.add(phase0Field, gbc);

            gbc.gridx = 3; filterPairsPanel.add(createLabel("Mag 0:"), gbc);
            JTextField mag0Field = createTextField(String.valueOf(currentlyEditedTone.filter.magnitudes[1][0][i]));
            mag0Field.getDocument().addDocumentListener(new DocumentListener() {
                public void changedUpdate(DocumentEvent e) { updateFilterMagnitude(1, 0, pairIndex, mag0Field); }
                public void removeUpdate(DocumentEvent e) { updateFilterMagnitude(1, 0, pairIndex, mag0Field); }
                public void insertUpdate(DocumentEvent e) { updateFilterMagnitude(1, 0, pairIndex, mag0Field); }
            });
            gbc.gridx = 4; filterPairsPanel.add(mag0Field, gbc);

            gbc.gridx = 1; gbc.gridy++; filterPairsPanel.add(createLabel("Phase 1:"), gbc);
            JTextField phase1Field = createTextField(String.valueOf(currentlyEditedTone.filter.phases[1][1][i]));
            phase1Field.getDocument().addDocumentListener(new DocumentListener() {
                public void changedUpdate(DocumentEvent e) { updateFilterPhase(1, 1, pairIndex, phase1Field); }
                public void removeUpdate(DocumentEvent e) { updateFilterPhase(1, 1, pairIndex, phase1Field); }
                public void insertUpdate(DocumentEvent e) { updateFilterPhase(1, 1, pairIndex, phase1Field); }
            });
            gbc.gridx = 2; filterPairsPanel.add(phase1Field, gbc);

            gbc.gridx = 3; filterPairsPanel.add(createLabel("Mag 1:"), gbc);
            JTextField mag1Field = createTextField(String.valueOf(currentlyEditedTone.filter.magnitudes[1][1][i]));
            mag1Field.getDocument().addDocumentListener(new DocumentListener() {
                public void changedUpdate(DocumentEvent e) { updateFilterMagnitude(1, 1, pairIndex, mag1Field); }
                public void removeUpdate(DocumentEvent e) { updateFilterMagnitude(1, 1, pairIndex, mag1Field); }
                public void insertUpdate(DocumentEvent e) { updateFilterMagnitude(1, 1, pairIndex, mag1Field); }
            });
            gbc.gridx = 4; filterPairsPanel.add(mag1Field, gbc);
        }

        filterPairsPanel.revalidate();
        filterPairsPanel.repaint();
    }

    private void updateFilterPhase(int direction, int phaseType, int pairIndex, JTextField field) {
        if (currentlyEditedTone != null && currentlyEditedTone.filter != null) {
            try {
                // Ensure index is within bounds before assignment
                if (pairIndex < currentlyEditedTone.filter.phases[direction][phaseType].length) {
                    currentlyEditedTone.filter.phases[direction][phaseType][pairIndex] = Integer.parseInt(field.getText());
                    debounceTimer.restart();
                }
            } catch (NumberFormatException e) {
                throw new RuntimeException(e);
            }
        }
    }

    private void updateFilterMagnitude(int direction, int phaseType, int pairIndex, JTextField field) {
        if (currentlyEditedTone != null && currentlyEditedTone.filter != null) {
            try {
                // Ensure index is within bounds before assignment
                if (pairIndex < currentlyEditedTone.filter.magnitudes[direction][phaseType].length) {
                    currentlyEditedTone.filter.magnitudes[direction][phaseType][pairIndex] = Integer.parseInt(field.getText());
                    debounceTimer.restart();
                }
            } catch (NumberFormatException e) {
                throw new RuntimeException(e);
            }
        }
    }


    // Add all listeners once during initialization
    private void addGlobalParameterListeners() {
        // ChangeListener for sliders - triggers debounce timer
        ChangeListener sliderChangeListener = _ -> {
            if (currentlyEditedTone != null) {
                currentlyEditedTone.duration = durationSlider.getValue();
                currentlyEditedTone.offset = offsetSlider.getValue();
                currentlyEditedTone.delayTime = delayTimeSlider.getValue();
                currentlyEditedTone.delayDecay = delayDecaySlider.getValue();

                // Update selected oscillator parameters
                if (currentlySelectedOscillatorIndex >= 0 && currentlySelectedOscillatorIndex < MAX_OSCILLATORS) {
                    currentlyEditedTone.oscillatorVolume[currentlySelectedOscillatorIndex] = oscVolumeSlider.getValue();
                    currentlyEditedTone.oscillatorPitch[currentlySelectedOscillatorIndex] = oscPitchSlider.getValue();
                    currentlyEditedTone.oscillatorDelays[currentlySelectedOscillatorIndex] = oscDelaySlider.getValue();
                }

                if (currentlyEditedTone.filter != null) {
                    currentlyEditedTone.filter.unity[0] = filterUnity0Slider.getValue();
                    currentlyEditedTone.filter.unity[1] = filterUnity1Slider.getValue();
                    // Update the filter's internal pairs directly from the slider values
                    currentlyEditedTone.filter.pairs[0] = filterPairs0Slider.getValue();
                    currentlyEditedTone.filter.pairs[1] = filterPairs1Slider.getValue();
                    updateFilterPairsPanel(filterEnabledCheckbox.isSelected()); // Rebuild filter pair UI
                }
                debounceTimer.restart(); // Restart the timer on every change
            }
        };

        durationSlider.addChangeListener(sliderChangeListener);
        offsetSlider.addChangeListener(sliderChangeListener);
        delayTimeSlider.addChangeListener(sliderChangeListener);
        delayDecaySlider.addChangeListener(sliderChangeListener);
        oscVolumeSlider.addChangeListener(sliderChangeListener);
        oscPitchSlider.addChangeListener(sliderChangeListener);
        oscDelaySlider.addChangeListener(sliderChangeListener);
        filterUnity0Slider.addChangeListener(sliderChangeListener);
        filterUnity1Slider.addChangeListener(sliderChangeListener);
        filterPairs0Slider.addChangeListener(sliderChangeListener);
        filterPairs1Slider.addChangeListener(sliderChangeListener);


        // Helper to add listeners to an envelope's UI components
        java.util.function.BiConsumer<SoundEnvelope, JComboBox<String>> formListener = (env, combo) -> combo.addActionListener(_ -> {
            if (currentlyEditedTone != null && env != null) { // Ensure envelope object exists
                int selectedIndex = combo.getSelectedIndex();
                env.form = selectedIndex == 4 ? 0 : selectedIndex + 1; // "Off" is 0, others 1-4
                debounceTimer.restart();
            }
        });
        java.util.function.Consumer<JTextField> createTextFieldListener = (field) -> field.getDocument().addDocumentListener(new DocumentListener() {
            private void updateValue() {
                if (currentlyEditedTone != null) {
                    try {
                        int parsedValue;
                        String text = field.getText();
                        if (text.isEmpty()) {
                            parsedValue = 0;
                        } else {
                            parsedValue = Integer.parseInt(text);
                        }

                        // Update the specific field based on its reference
                        if (field == pitchEnvelopeStartField) currentlyEditedTone.pitch.start = parsedValue;
                        else if (field == pitchEnvelopeEndField) currentlyEditedTone.pitch.end = parsedValue;
                        else if (field == volumeEnvelopeStartField) currentlyEditedTone.volume.start = parsedValue;
                        else if (field == volumeEnvelopeEndField) currentlyEditedTone.volume.end = parsedValue;
                            // For optional envelopes, check if the object itself is not null before accessing its fields
                        else if (field == pitchModifierStartField && currentlyEditedTone.pitchModifier != null) currentlyEditedTone.pitchModifier.start = parsedValue;
                        else if (field == pitchModifierEndField && currentlyEditedTone.pitchModifier != null) currentlyEditedTone.pitchModifier.end = parsedValue;
                        else if (field == pitchModifierAmplitudeStartField && currentlyEditedTone.pitchModifierAmplitude != null) currentlyEditedTone.pitchModifierAmplitude.start = parsedValue;
                        else if (field == pitchModifierAmplitudeEndField && currentlyEditedTone.pitchModifierAmplitude != null) currentlyEditedTone.pitchModifierAmplitude.end = parsedValue;
                        else if (field == volumeMultiplierStartField && currentlyEditedTone.volumeMultiplier != null) currentlyEditedTone.volumeMultiplier.start = parsedValue;
                        else if (field == volumeMultiplierEndField && currentlyEditedTone.volumeMultiplier != null) currentlyEditedTone.volumeMultiplier.end = parsedValue;
                        else if (field == volumeMultiplierAmplitudeStartField && currentlyEditedTone.volumeMultiplierAmplitude != null) currentlyEditedTone.volumeMultiplierAmplitude.start = parsedValue;
                        else if (field == volumeMultiplierAmplitudeEndField && currentlyEditedTone.volumeMultiplierAmplitude != null) currentlyEditedTone.volumeMultiplierAmplitude.end = parsedValue;
                        else if (field == releaseStartField && currentlyEditedTone.release != null) currentlyEditedTone.release.start = parsedValue;
                        else if (field == releaseEndField && currentlyEditedTone.release != null) currentlyEditedTone.release.end = parsedValue;
                        else if (field == attackStartField && currentlyEditedTone.attack != null) currentlyEditedTone.attack.start = parsedValue;
                        else if (field == attackEndField && currentlyEditedTone.attack != null) currentlyEditedTone.attack.end = parsedValue;

                        // Only restart timer if parsing was successful
                        debounceTimer.restart();

                        // Optional: Reset field border to default if it was previously red
                        field.setBorder(UIManager.getBorder("TextField.border")); // Or your default border
                    } catch (NumberFormatException ex) {
                        // Provide immediate visual feedback to the user
                        // 1. Change border to red
                        field.setBorder(BorderFactory.createLineBorder(Color.RED, 1));
                        // 2. Show a small tooltip or status message (optional)
                        // field.setToolTipText("Please enter a valid number.");
                        // 3. You can also print to console for debugging, but user feedback is key
                        System.out.println("Invalid number: " + field.getText() + " - " + ex.getMessage());

                        // IMPORTANT: Do NOT restart debounceTimer here, as the input is invalid
                    }
                }
            }
            @Override public void insertUpdate(DocumentEvent e) { updateValue(); }
            @Override public void removeUpdate(DocumentEvent e) { updateValue(); }
            @Override public void changedUpdate(DocumentEvent e) { updateValue(); }
        });

        // Apply listeners to all envelope text fields
        createTextFieldListener.accept(pitchEnvelopeStartField);
        createTextFieldListener.accept(pitchEnvelopeEndField);
        createTextFieldListener.accept(volumeEnvelopeStartField);
        createTextFieldListener.accept(volumeEnvelopeEndField);
        createTextFieldListener.accept(pitchModifierStartField);
        createTextFieldListener.accept(pitchModifierEndField);
        createTextFieldListener.accept(pitchModifierAmplitudeStartField);
        createTextFieldListener.accept(pitchModifierAmplitudeEndField);
        createTextFieldListener.accept(volumeMultiplierStartField);
        createTextFieldListener.accept(volumeMultiplierEndField);
        createTextFieldListener.accept(volumeMultiplierAmplitudeStartField);
        createTextFieldListener.accept(volumeMultiplierAmplitudeEndField);
        createTextFieldListener.accept(releaseStartField);
        createTextFieldListener.accept(releaseEndField);
        createTextFieldListener.accept(attackStartField);
        createTextFieldListener.accept(attackEndField);


        // Apply form listeners to all envelopes
        // For mandatory envelopes, they are always present, so we can directly pass their instances
        formListener.accept(currentlyEditedTone.pitch, pitchEnvelopeFormComboBox);
        formListener.accept(currentlyEditedTone.volume, volumeEnvelopeFormComboBox);
        // For optional envelopes, add listeners only if they are not null (i.e., created by enabling checkbox)
        // The checkbox action listener will handle creating the object if it's null.
        pitchModifierFormComboBox.addActionListener(_ -> {
            if (currentlyEditedTone != null && currentlyEditedTone.pitchModifier != null) {
                int selectedIndex = pitchModifierFormComboBox.getSelectedIndex();
                currentlyEditedTone.pitchModifier.form = selectedIndex == 4 ? 0 : selectedIndex + 1;
                debounceTimer.restart();
            }
        });
        pitchModifierAmplitudeFormComboBox.addActionListener(_ -> {
            if (currentlyEditedTone != null && currentlyEditedTone.pitchModifierAmplitude != null) {
                int selectedIndex = pitchModifierAmplitudeFormComboBox.getSelectedIndex();
                currentlyEditedTone.pitchModifierAmplitude.form = selectedIndex == 4 ? 0 : selectedIndex + 1;
                debounceTimer.restart();
            }
        });
        volumeMultiplierFormComboBox.addActionListener(_ -> {
            if (currentlyEditedTone != null && currentlyEditedTone.volumeMultiplier != null) {
                int selectedIndex = volumeMultiplierFormComboBox.getSelectedIndex();
                currentlyEditedTone.volumeMultiplier.form = selectedIndex == 4 ? 0 : selectedIndex + 1;
                debounceTimer.restart();
            }
        });
        volumeMultiplierAmplitudeFormComboBox.addActionListener(_ -> {
            if (currentlyEditedTone != null && currentlyEditedTone.volumeMultiplierAmplitude != null) {
                int selectedIndex = volumeMultiplierAmplitudeFormComboBox.getSelectedIndex();
                currentlyEditedTone.volumeMultiplierAmplitude.form = selectedIndex == 4 ? 0 : selectedIndex + 1;
                debounceTimer.restart();
            }
        });
        releaseFormComboBox.addActionListener(_ -> {
            if (currentlyEditedTone != null && currentlyEditedTone.release != null) {
                int selectedIndex = releaseFormComboBox.getSelectedIndex();
                currentlyEditedTone.release.form = selectedIndex == 4 ? 0 : selectedIndex + 1;
                debounceTimer.restart();
            }
        });
        attackFormComboBox.addActionListener(_ -> {
            if (currentlyEditedTone != null && currentlyEditedTone.attack != null) {
                int selectedIndex = attackFormComboBox.getSelectedIndex();
                currentlyEditedTone.attack.form = selectedIndex == 4 ? 0 : selectedIndex + 1;
                debounceTimer.restart();
            }
        });

        pitchModifierEnabledCheckbox.addActionListener(_ -> {
            if (currentlyEditedTone != null) {
                boolean enabled = pitchModifierEnabledCheckbox.isSelected();
                if (enabled) {
                    if (currentlyEditedTone.pitchModifier == null) {
                        currentlyEditedTone.pitchModifier = new SoundEnvelope();
                        currentlyEditedTone.pitchModifierAmplitude = new SoundEnvelope();
                        currentlyEditedTone.pitchModifier.form = 1;
                        currentlyEditedTone.pitchModifierAmplitude.form = 1;
                    }
                } else {
                    currentlyEditedTone.pitchModifier = null;
                    currentlyEditedTone.pitchModifierAmplitude = null;
                }
                updateEnvelopeUI(currentlyEditedTone.pitchModifier, pitchModifierFormComboBox, pitchModifierStartField, pitchModifierEndField, pitchModifierEnvelopeGraphPanel, pitchModifierEnabledCheckbox);
                updateEnvelopeUI(currentlyEditedTone.pitchModifierAmplitude, pitchModifierAmplitudeFormComboBox, pitchModifierAmplitudeStartField, pitchModifierAmplitudeEndField, pitchModifierAmplitudeEnvelopeGraphPanel, pitchModifierAmplitudeEnabledCheckbox);
                debounceTimer.restart();
            }
        });

        pitchModifierAmplitudeEnabledCheckbox.addActionListener(_ -> {
            if (currentlyEditedTone != null) {
                boolean enabled = pitchModifierAmplitudeEnabledCheckbox.isSelected();
                if (enabled) {
                    if (currentlyEditedTone.pitchModifierAmplitude == null) {
                        currentlyEditedTone.pitchModifierAmplitude = new SoundEnvelope();
                        currentlyEditedTone.pitchModifierAmplitude.form = 1;
                    }
                } else {
                    currentlyEditedTone.pitchModifierAmplitude = null;
                }
                updateEnvelopeUI(currentlyEditedTone.pitchModifierAmplitude, pitchModifierAmplitudeFormComboBox, pitchModifierAmplitudeStartField, pitchModifierAmplitudeEndField, pitchModifierAmplitudeEnvelopeGraphPanel, pitchModifierAmplitudeEnabledCheckbox);
                debounceTimer.restart();
            }
        });

        volumeMultiplierEnabledCheckbox.addActionListener(_ -> {
            if (currentlyEditedTone != null) {
                boolean enabled = volumeMultiplierEnabledCheckbox.isSelected();
                if (enabled) {
                    if (currentlyEditedTone.volumeMultiplier == null) {
                        currentlyEditedTone.volumeMultiplier = new SoundEnvelope();
                        currentlyEditedTone.volumeMultiplierAmplitude = new SoundEnvelope();
                        currentlyEditedTone.volumeMultiplier.form = 1;
                        currentlyEditedTone.volumeMultiplierAmplitude.form = 1;
                    }
                } else {
                    currentlyEditedTone.volumeMultiplier = null;
                    currentlyEditedTone.volumeMultiplierAmplitude = null;
                }
                updateEnvelopeUI(currentlyEditedTone.volumeMultiplier, volumeMultiplierFormComboBox, volumeMultiplierStartField, volumeMultiplierEndField, volumeMultiplierEnvelopeGraphPanel, volumeMultiplierEnabledCheckbox);
                updateEnvelopeUI(currentlyEditedTone.volumeMultiplierAmplitude, volumeMultiplierAmplitudeFormComboBox, volumeMultiplierAmplitudeStartField, volumeMultiplierAmplitudeEndField, volumeMultiplierAmplitudeEnvelopeGraphPanel, volumeMultiplierAmplitudeEnabledCheckbox);
                debounceTimer.restart();
            }
        });

        volumeMultiplierAmplitudeEnabledCheckbox.addActionListener(_ -> {
            if (currentlyEditedTone != null) {
                boolean enabled = volumeMultiplierAmplitudeEnabledCheckbox.isSelected();
                if (enabled) {
                    if (currentlyEditedTone.volumeMultiplierAmplitude == null) {
                        currentlyEditedTone.volumeMultiplierAmplitude = new SoundEnvelope();
                        currentlyEditedTone.volumeMultiplierAmplitude.form = 1;
                    }
                } else {
                    currentlyEditedTone.volumeMultiplierAmplitude = null;
                }
                updateEnvelopeUI(currentlyEditedTone.volumeMultiplierAmplitude, volumeMultiplierAmplitudeFormComboBox, volumeMultiplierAmplitudeStartField, volumeMultiplierAmplitudeEndField, volumeMultiplierAmplitudeEnvelopeGraphPanel, volumeMultiplierAmplitudeEnabledCheckbox);
                debounceTimer.restart();
            }
        });

        releaseEnabledCheckbox.addActionListener(_ -> {
            if (currentlyEditedTone != null) {
                boolean enabled = releaseEnabledCheckbox.isSelected();
                if (enabled) {
                    if (currentlyEditedTone.release == null) {
                        currentlyEditedTone.release = new SoundEnvelope();
                        currentlyEditedTone.attack = new SoundEnvelope(); // Attack is paired with Release
                        currentlyEditedTone.release.form = 1;
                        currentlyEditedTone.attack.form = 1;
                    }
                } else {
                    currentlyEditedTone.release = null;
                    currentlyEditedTone.attack = null;
                }
                updateEnvelopeUI(currentlyEditedTone.release, releaseFormComboBox, releaseStartField, releaseEndField, releaseEnvelopeGraphPanel, releaseEnabledCheckbox);
                updateEnvelopeUI(currentlyEditedTone.attack, attackFormComboBox, attackStartField, attackEndField, attackEnvelopeGraphPanel, attackEnabledCheckbox);
                debounceTimer.restart();
            }
        });
        attackEnabledCheckbox.addActionListener(_ -> {
            if (currentlyEditedTone != null) {
                boolean enabled = attackEnabledCheckbox.isSelected();
                if (enabled && currentlyEditedTone.release != null) {
                    if (currentlyEditedTone.attack == null) {
                        currentlyEditedTone.attack = new SoundEnvelope();
                        currentlyEditedTone.attack.form = 1;
                    }
                } else {
                    currentlyEditedTone.attack = null;
                }
                updateEnvelopeUI(currentlyEditedTone.attack, attackFormComboBox, attackStartField, attackEndField, attackEnvelopeGraphPanel, attackEnabledCheckbox);
                debounceTimer.restart();
            }
        });


        filterEnabledCheckbox.addActionListener(_ -> {
            if (currentlyEditedTone != null) {
                boolean enabled = filterEnabledCheckbox.isSelected();
                if (enabled) {
                    if (currentlyEditedTone.filter == null) {
                        currentlyEditedTone.filter = new SoundFilter();
                        currentlyEditedTone.filterEnvelope = new SoundEnvelope();
                    }
                    currentlyEditedTone.filter.pairs[0] = 1;
                    currentlyEditedTone.filter.pairs[1] = 1;
                    currentlyEditedTone.filter.unity[0] = 32768;
                    currentlyEditedTone.filter.unity[1] = 32768;
                } else {
                    currentlyEditedTone.filter = null;
                    currentlyEditedTone.filterEnvelope = null;
                }

                filterUnity0Slider.setEnabled(enabled);
                filterUnity1Slider.setEnabled(enabled);
                filterPairs0Slider.setEnabled(enabled);
                filterPairs1Slider.setEnabled(enabled);

                if (currentlyEditedTone.filter != null) {
                    filterUnity0Slider.setValue(currentlyEditedTone.filter.unity[0]);
                    filterUnity1Slider.setValue(currentlyEditedTone.filter.unity[1]);
                    filterPairs0Slider.setValue(currentlyEditedTone.filter.pairs[0]);
                    filterPairs1Slider.setValue(currentlyEditedTone.filter.pairs[1]);
                } else {
                    filterUnity0Slider.setValue(0);
                    filterUnity1Slider.setValue(0);
                    filterPairs0Slider.setValue(0);
                    filterPairs1Slider.setValue(0);
                }

                updateFilterPairsPanel(enabled); // Rebuild filter pair UI
                debounceTimer.restart();
            }
        });
    }

    private void updateWaveformDisplay() {
        if (currentlyEditedTone != null) {
            try {
                // Ensure duration is positive to avoid issues with zero-length arrays
                int durationMs = currentlyEditedTone.duration > 0 ? currentlyEditedTone.duration : 1;
                int steps = durationMs * 22050 / 1000;

                if (steps <= 0) {
                    waveformDisplayPanel.setAudioData(new byte[0]);
                    return;
                }
                // Synthesize returns int[] (16-bit samples)
                int[] rawSamplesInt = currentlyEditedTone.synthesize(steps, durationMs);

                // Convert to byte array based on selected output bit depth for display
                // Note: The mix() method in SoundEffect converts to 8-bit.
                // For individual tone display, we convert the 16-bit int[] to the selected bit depth.
                byte[] audioDataForDisplay = convertRawSamplesToPlaybackFormat(rawSamplesInt, outputBitDepth);
                waveformDisplayPanel.setAudioData(audioDataForDisplay);
            } catch (Exception ex) {
                System.out.println("Error synthesizing tone for waveform display: " + ex.getMessage());
                waveformDisplayPanel.setAudioData(new byte[0]);
            }
        } else {
            waveformDisplayPanel.setAudioData(new byte[0]);
        }
    }

    /**
     * Converts an array of 16-bit integer samples to a byte array
     * with the specified bit depth (8 or 16) for playback or display.
     * The input `samples` array is assumed to contain 16-bit signed values
     * (e.g., -32768 to 32767) stored in Java `int` primitives.
     *
     * @param samples The input array of integer samples (16-bit).
     * @param targetBitDepth The desired output bit depth (8 or 16).
     * @return A new byte array containing the converted audio data.
     */
    private byte[] convertRawSamplesToPlaybackFormat(int[] samples, int targetBitDepth) {
        if (targetBitDepth == 8) {
            byte[] audioData = new byte[samples.length];
            for (int i = 0; i < samples.length; ++i) {
                // Convert 16-bit sample to 8-bit by taking the high byte
                int sample = samples[i] >> 8;
                // Clamp to -128 to 127 range for 8-bit signed byte
                if ((sample + 128 & -256) != 0) { // This is the original clamping logic from SoundEffect.mix
                    sample = sample >> 31 ^ 127; // Equivalent to (sample < -128) ? -128 : ((sample > 127) ? 127 : sample)
                }
                audioData[i] = (byte) sample;
            }
            return audioData;
        } else if (targetBitDepth == 16) {
            byte[] audioData = new byte[samples.length * 2]; // 2 bytes per 16-bit sample
            for (int i = 0; i < samples.length; ++i) {
                int sample = samples[i]; // Sample is already in 16-bit range (-32768 to 32767)
                // Write as little-endian: LSB first, then MSB
                audioData[i * 2] = (byte) (sample & 0xFF);         // Low byte
                audioData[i * 2 + 1] = (byte) ((sample >> 8) & 0xFF); // High byte
            }
            return audioData;
        } else {
            throw new IllegalArgumentException("Unsupported bit depth: " + targetBitDepth + ". Only 8 or 16 are supported.");
        }
    }


    // --- Sound Playback ---

    private void playSound(byte[] audioData, int sampleRate, int bitDepth) {
        stopSoundEffect(); // Stop any currently playing sound

        if (audioData.length == 0) {
            JOptionPane.showMessageDialog(this, "Generated audio data is empty.", "Playback Error", JOptionPane.WARNING_MESSAGE);
            return;
        }

        try {
            // AudioFormat: sampleRate, bitDepth, mono (1 channel), signed, little-endian
            AudioFormat format = new AudioFormat(sampleRate, bitDepth, 1, true, false);

            DataLine.Info info = new DataLine.Info(SourceDataLine.class, format);
            if (!AudioSystem.isLineSupported(info)) {
                JOptionPane.showMessageDialog(this, "Audio line not supported for format: " + format, "Playback Error", JOptionPane.ERROR_MESSAGE);
                return;
            }

            SourceDataLine line = (SourceDataLine) AudioSystem.getLine(info);
            line.open(format);
            line.start();

            // Use a separate thread for playback to avoid freezing the UI
            Thread playbackThread = new Thread(() -> {
                try {
                    // Create a Clip to manage playback (allows stopping)
                    audioClip = AudioSystem.getClip();
                    audioClip.open(format, audioData, 0, audioData.length);
                    audioClip.start();
                    audioClip.addLineListener(event -> {
                        if (event.getType() == LineEvent.Type.STOP) {
                            audioClip.close();
                        }
                    });

                } catch (LineUnavailableException e) {
                    JOptionPane.showMessageDialog(this, "Error playing sound: " + e.getMessage(), "Playback Error", JOptionPane.ERROR_MESSAGE);
                }
            });
            playbackThread.start();

        } catch (Exception ex) {
            JOptionPane.showMessageDialog(this, "Error playing sound: " + ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
        }
    }

    private SoundEffect newSoundEffectInstance() {
        SoundEffect newSoundEffect = new SoundEffect();

        for (int i = 0; i < currentSoundTones.size(); i++) {
            if (i < newSoundEffect.soundTones.length) {
                newSoundEffect.soundTones[i] = currentSoundTones.get(i);
            }
        }

        for (int i = currentSoundTones.size(); i < newSoundEffect.soundTones.length; i++) {
            newSoundEffect.soundTones[i] = null;
        }

        int maxDuration = 0;
        for (SoundTone tone : currentSoundTones) {
            if (tone != null) {
                maxDuration = Math.max(maxDuration, tone.duration + tone.offset);
            }
        }
        newSoundEffect.start = 0;
        newSoundEffect.end = maxDuration;
        return newSoundEffect;
    }

    private void playSoundEffect() {
        if (currentSoundTones.isEmpty()) {
            JOptionPane.showMessageDialog(this, "No Sound Tones to play.", "Playback Error", JOptionPane.WARNING_MESSAGE);
            return;
        }

        try {
            SoundEffect soundEffectToPlay = newSoundEffectInstance();
            byte[] rawSamples8Bit = soundEffectToPlay.toRawSound().audioDataByte;
            byte[] audioDataForPlayback;
            if (outputBitDepth == 8) {
                audioDataForPlayback = rawSamples8Bit;
            } else {
                audioDataForPlayback = new byte[rawSamples8Bit.length * 2];
                for (int i = 0; i < rawSamples8Bit.length; i++) {
                    short sample16 = (short) (rawSamples8Bit[i] * 256);
                    audioDataForPlayback[i * 2] = (byte) (sample16 & 0xFF);
                    audioDataForPlayback[i * 2 + 1] = (byte) ((sample16 >> 8) & 0xFF);
                }
            }

            int sampleRate = 22050;

            playSound(audioDataForPlayback, sampleRate, outputBitDepth);

        } catch (Exception ex) {
            JOptionPane.showMessageDialog(this, "Error generating or playing full sound effect: " + ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
        }
    }

    private void playSelectedSoundTone() {
        int selectedIndex = soundToneList.getSelectedIndex();
        if (selectedIndex == -1) {
            JOptionPane.showMessageDialog(this, "Please select a Sound Tone to play.", "No Selection", JOptionPane.WARNING_MESSAGE);
            return;
        }

        SoundTone selectedTone = currentSoundTones.get(selectedIndex);

        try {
            int steps = selectedTone.duration * 22050 / 1000;
            if (steps <= 0) {
                JOptionPane.showMessageDialog(this, "Selected tone has zero or negative duration, cannot synthesize.", "Synthesis Error", JOptionPane.WARNING_MESSAGE);
                return;
            }
            int[] rawSamplesInt = selectedTone.synthesize(steps, selectedTone.duration);
            byte[] audioData = convertRawSamplesToPlaybackFormat(rawSamplesInt, outputBitDepth);
            playSound(audioData, 22050, outputBitDepth);
        } catch (Exception ex) {
            JOptionPane.showMessageDialog(this, "Error generating or playing selected tone: " + ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
        }
    }


    private void stopSoundEffect() {
        if (audioClip != null && audioClip.isRunning()) {
            audioClip.stop();
            audioClip.close();
        }
    }

    /**
     * Handles saving the current sound tones to a file in a format
     * that would theoretically be compliant with the JagFX SoundEffect decoder.
     */
    private void saveJagFXFile() {
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setDialogTitle("Save Sound Effect as JagFX File");
        fileChooser.setFileFilter(new FileNameExtensionFilter("JagFX Sound Effect Files (*.jfx)", "jfx"));
        fileChooser.setSelectedFile(new File("my_sound_effect.jfx"));

        int userSelection = fileChooser.showSaveDialog(this);

        if (userSelection == JFileChooser.APPROVE_OPTION) {
            File fileToSave = fileChooser.getSelectedFile();
            // Ensure the file has the correct extension
            if (!fileToSave.getName().toLowerCase().endsWith(".jfx")) {
                fileToSave = new File(fileToSave.getAbsolutePath() + ".jfx");
            }

            try (FileOutputStream fos = new FileOutputStream(fileToSave);
                 DataOutputStream dos = new DataOutputStream(fos)) {

                // Encode the SoundEffect into a byte array
                byte[] encodedData = newSoundEffectInstance().encode();

                // Write the collected bytes to the file
                dos.write(encodedData);

                System.out.println("Sound effect saved to: " + fileToSave.getAbsolutePath());
                JOptionPane.showMessageDialog(this, "Sound effect saved successfully!", "Save Complete", JOptionPane.INFORMATION_MESSAGE);

            } catch (IOException ex) {
                JOptionPane.showMessageDialog(this, "Error saving file: " + ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
            }
        }
    }

    /**
     * Handles loading a sound effect from a JagFX file.
     */
    private void importJagFXFile() {
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setDialogTitle("Import Sound Effect from JagFX File");
        fileChooser.setFileFilter(new FileNameExtensionFilter("JagFX Sound Effect Files (*.jfx), (*.dat)", "jfx", "dat"));

        int userSelection = fileChooser.showOpenDialog(this);

        if (userSelection == JFileChooser.APPROVE_OPTION) {
            File fileToLoad = fileChooser.getSelectedFile();

            try (FileInputStream fis = new FileInputStream(fileToLoad);
                 BufferedInputStream bis = new BufferedInputStream(fis)) {

                byte[] fileBytes = bis.readAllBytes();
                Buffer inputBuffer = new Buffer(fileBytes);
                SoundEffect loadedSoundEffect = new SoundEffect(inputBuffer);

                currentSoundTones.clear();
                for (SoundTone tone : loadedSoundEffect.soundTones) {
                    if (tone != null) {
                        currentSoundTones.add(tone);
                    }
                }

                if (currentSoundTones.isEmpty()) {
                    currentSoundTones.add(new SoundTone());
                    JOptionPane.showMessageDialog(this, "No valid sound tones found in the file. A default tone has been loaded.", "Import Warning", JOptionPane.WARNING_MESSAGE);
                }

                updateToneList();
                JOptionPane.showMessageDialog(this, "Sound effect loaded successfully!", "Import Complete", JOptionPane.INFORMATION_MESSAGE);

            } catch (FileNotFoundException ex) {
                JOptionPane.showMessageDialog(this, "File not found: " + ex.getMessage(), "Import Error", JOptionPane.ERROR_MESSAGE);
            } catch (IOException ex) {
                JOptionPane.showMessageDialog(this, "Error reading file: " + ex.getMessage(), "Import Error", JOptionPane.ERROR_MESSAGE);
            } catch (ArrayIndexOutOfBoundsException ex) {
                JOptionPane.showMessageDialog(this, "Error decoding file: The file format appears to be corrupted or incompatible. " + "Details: " + ex.getMessage() + "\n\nPlease ensure it's a valid JagFX sound effect file.", "Import Error: File Format Mismatch", JOptionPane.ERROR_MESSAGE);
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(this, "An unexpected error occurred during import: " + ex.getMessage(), "Import Error", JOptionPane.ERROR_MESSAGE);
            }
        }
    }

    /**
     * Exports the current full sound effect to a WAV file.
     */
    private void exportToWavFile() {
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setDialogTitle("Export Sound Effect to WAV File");
        fileChooser.setFileFilter(new FileNameExtensionFilter("WAV Audio Files (*.wav)", "wav"));
        fileChooser.setSelectedFile(new File("exported_sound_effect.wav"));

        int userSelection = fileChooser.showSaveDialog(this);

        if (userSelection == JFileChooser.APPROVE_OPTION) {
            File fileToSave = fileChooser.getSelectedFile();
            if (!fileToSave.getName().toLowerCase().endsWith(".wav")) {
                fileToSave = new File(fileToSave.getAbsolutePath() + ".wav");
            }

            try {
                SoundEffect tempSoundEffect = newSoundEffectInstance();
                byte[] combinedRawSamples8Bit = tempSoundEffect.toRawSound().audioDataByte;
                byte[] wavData16Bit = new byte[combinedRawSamples8Bit.length * 2];
                for (int i = 0; i < combinedRawSamples8Bit.length; i++) {
                    short sample16 = (short) (combinedRawSamples8Bit[i] * 256);
                    wavData16Bit[i * 2] = (byte) (sample16 & 0xFF);
                    wavData16Bit[i * 2 + 1] = (byte) ((sample16 >> 8) & 0xFF);
                }

                AudioFormat format = new AudioFormat(22050, 16, 1, true, false);

                ByteArrayInputStream bais = new ByteArrayInputStream(wavData16Bit);
                AudioInputStream audioInputStream = new AudioInputStream(bais, format, wavData16Bit.length / format.getFrameSize());

                AudioSystem.write(audioInputStream, AudioFileFormat.Type.WAVE, fileToSave);

                audioInputStream.close();
                bais.close();

                JOptionPane.showMessageDialog(this, "Sound effect exported to WAV successfully!", "Export Complete", JOptionPane.INFORMATION_MESSAGE);

            } catch (IOException ex) {
                JOptionPane.showMessageDialog(this, "Error exporting WAV file: " + ex.getMessage(), "Export Error", JOptionPane.ERROR_MESSAGE);
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(this, "An unexpected error occurred during WAV export: " + ex.getMessage(), "Export Error", JOptionPane.ERROR_MESSAGE);
            }
        }
    }


    public static void main(String[] args) {
        SwingUtilities.invokeLater(SoundEffectEditor::new);
    }

    private class EnvelopeGraphPanel extends JPanel {

        private SoundEnvelope envelope;
        private final Color gridColor;
        private final Color lineColor;

        public EnvelopeGraphPanel(String title, Color bgColor, Color gridColor, Color lineColor) {
            this.gridColor = gridColor;
            this.lineColor = lineColor;
            setBackground(bgColor);
            setBorder(BorderFactory.createTitledBorder(
                    BorderFactory.createLineBorder(BORDER_COLOR),
                    title,
                    TitledBorder.LEFT,
                    TitledBorder.TOP,
                    new Font("Arial", Font.BOLD, 12),
                    HIGHLIGHT_COLOR_GREEN // Changed to green highlight
            ));
        }

        public void setEnvelope(SoundEnvelope envelope) {
            this.envelope = envelope;
            repaint();
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g;
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            int width = getWidth();
            int height = getHeight();

            // Draw grid
            g2d.setColor(gridColor);
            for (int i = 0; i <= 10; i++) {
                g2d.drawLine(0, i * height / 10, width, i * height / 10); // Horizontal
                g2d.drawLine(i * width / 10, 0, i * width / 10, height); // Vertical
            }

            if (envelope != null && envelope.form != 0) { // Only draw if envelope is active
                g2d.setColor(lineColor);
                g2d.setStroke(new BasicStroke(2)); // Thicker line for emphasis
                Path2D.Float path = new Path2D.Float();

                // To draw a more accurate envelope, we need to simulate its doStep method
                // over a range of 'period' (time).
                // Let's simulate 1000 steps for the graph.
                int graphSteps = 1000;
                int maxEnvelopeValue = 32768; // Max value for envelope output (from doStep)

                // Create a temporary envelope instance to avoid modifying the actual one
                SoundEnvelope tempEnvelope = new SoundEnvelope();
                tempEnvelope.form = envelope.form;
                tempEnvelope.start = envelope.start;
                tempEnvelope.end = envelope.end;
                tempEnvelope.segments = envelope.segments; // Correctly copy the int value
                tempEnvelope.durations = Arrays.copyOf(envelope.durations, envelope.durations.length); // Correctly copy the array
                tempEnvelope.phases = Arrays.copyOf(envelope.phases, envelope.phases.length); // Correctly copy the array
                tempEnvelope.reset();

                // Calculate the period based on the current tone's duration
                int period = currentlyEditedTone != null && currentlyEditedTone.duration > 0 ? currentlyEditedTone.duration : 500; // Default if no tone selected or duration is 0

                for (int i = 0; i <= graphSteps; i++) {
                    float x = (float) i / graphSteps * width;
                    // Simulate the envelope's value at this point in time
                    int envelopeValue = tempEnvelope.doStep(period); // Use the tone's duration as period

                    // Normalize envelopeValue to 0-1 range (assuming 0-32768 output from doStep)
                    float normalizedValue = (float) envelopeValue / maxEnvelopeValue;
                    // Invert Y-axis for drawing (0 at top, max at bottom)
                    float y = height * (1.0f - normalizedValue);

                    if (i == 0) {
                        path.moveTo(x, y);
                    } else {
                        path.lineTo(x, y);
                    }
                }
                g2d.draw(path);
            }
        }
    }

    /**
     * Custom JPanel for drawing the waveform of audio data.
     */
    private class WaveformDisplayPanel extends JPanel {
        private byte[] audioData;
        private final Color waveformColor;

        public WaveformDisplayPanel(String title, Color bgColor, Color waveformColor) {
            this.waveformColor = waveformColor;
            setBackground(bgColor);
            setBorder(BorderFactory.createTitledBorder(
                    BorderFactory.createLineBorder(BORDER_COLOR),
                    title,
                    TitledBorder.LEFT,
                    TitledBorder.TOP,
                    new Font("Arial", Font.BOLD, 12),
                    HIGHLIGHT_COLOR_GREEN // Changed to green highlight
            ));
        }

        public void setAudioData(byte[] audioData) {
            this.audioData = audioData;
            repaint();
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g;
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            int width = getWidth();
            int height = getHeight();

            if (audioData != null && audioData.length > 0) {
                g2d.setColor(waveformColor);
                g2d.setStroke(new BasicStroke(1)); // Thin line for waveform
                Path2D.Float path = new Path2D.Float();

                int sampleSize = outputBitDepth / 8;
                int numSamples = audioData.length / sampleSize;

                if (numSamples == 0) return;

                float xScale = (float) width / numSamples;
                float yScale = height / 2.0f; // Center the waveform

                float firstSampleValue;
                if (outputBitDepth == 8) {
                    firstSampleValue = audioData[0];
                } else { // 16-bit
                    firstSampleValue = (short) (((audioData[1] & 0xFF) << 8) | (audioData[0] & 0xFF));
                }

                path.moveTo(0, yScale - (firstSampleValue * yScale / (outputBitDepth == 8 ? 128.0f : 32768.0f)));

                for (int i = 1; i < numSamples; i++) {
                    float x = i * xScale;
                    float sampleValue;

                    if (outputBitDepth == 8) {
                        sampleValue = audioData[i];
                    } else { // 16-bit
                        int byteIndex = i * 2;
                        sampleValue = (short) (((audioData[byteIndex + 1] & 0xFF) << 8) | (audioData[byteIndex] & 0xFF));
                    }

                    float y = yScale - (sampleValue * yScale / (outputBitDepth == 8 ? 128.0f : 32768.0f));
                    path.lineTo(x, y);
                }
                g2d.draw(path);
            }
        }
    }

    /**
     * Custom JPanel for drawing SoundFilter frequency response.
     * This is a highly simplified placeholder. A real filter graph
     * would require complex DSP calculations to determine frequency response.
     */
    private static class FilterGraphPanel extends JPanel {
        private SoundFilter filter;
        private final Color gridColor;
        private final Color lineColor;

        public FilterGraphPanel(String title, Color bgColor, Color gridColor, Color lineColor) {
            this.gridColor = gridColor;
            this.lineColor = lineColor;
            setBackground(bgColor);
            setBorder(BorderFactory.createTitledBorder(BorderFactory.createLineBorder(BORDER_COLOR), title, TitledBorder.LEFT, TitledBorder.TOP, new Font("Arial", Font.BOLD, 12), HIGHLIGHT_COLOR_GREEN));
        }

        public void setFilter(SoundFilter filter) {
            this.filter = filter;
            repaint();
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g;
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            int width = getWidth();
            int height = getHeight();

            // Draw grid
            g2d.setColor(gridColor);
            for (int i = 0; i <= 10; i++) {
                g2d.drawLine(0, i * height / 10, width, i * height / 10); // Horizontal
                g2d.drawLine(i * width / 10, 0, i * width / 10, height); // Vertical
            }

            if (filter != null && (filter.pairs[0] > 0 || filter.pairs[1] > 0)) {
                g2d.setColor(lineColor);
                g2d.setStroke(new BasicStroke(2)); // Thicker line for emphasis
                Path2D.Float path = new Path2D.Float();
                float unity0Normalized = (float) filter.unity[0] / 65535.0f;
                float unity1Normalized = (float) filter.unity[1] / 65535.0f;

                int startY = (int) (height * (1.0f - unity0Normalized));
                int endY = (int) (height * (1.0f - unity1Normalized));

                path.moveTo(0, startY);
                path.lineTo(width, endY);
                g2d.draw(path);
            }
        }
    }
}