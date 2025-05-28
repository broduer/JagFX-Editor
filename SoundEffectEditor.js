class BufferJS {
    constructor(data) {
        if (data instanceof Uint8Array || data instanceof ArrayBuffer) {
            this.arrayBuffer = (data instanceof ArrayBuffer) ? data : data.buffer;
            this.dataView = new DataView(this.arrayBuffer);
            this.readMode = true;
            this.length = this.arrayBuffer.byteLength;
        } else {
            // Write mode: initial capacity, can grow
            this.arrayBuffer = new ArrayBuffer(1024); // Initial capacity
            this.dataView = new DataView(this.arrayBuffer);
            this.readMode = false;
            this.length = 0; // Bytes written
        }
        this.offset = 0;
    }

    _ensureCapacity(needed) {
        if (this.offset + needed > this.arrayBuffer.byteLength) {
            const newCapacity = Math.max(this.arrayBuffer.byteLength * 2, this.offset + needed + 1024);
            const newBuffer = new ArrayBuffer(newCapacity);
            new Uint8Array(newBuffer).set(new Uint8Array(this.arrayBuffer, 0, this.length));
            this.arrayBuffer = newBuffer;
            this.dataView = new DataView(this.arrayBuffer);
        }
    }

    writeUnsignedByte(value) {
        if (this.readMode) throw new Error("Buffer is in read mode");
        this._ensureCapacity(1);
        this.dataView.setUint8(this.offset, value & 0xFF);
        this.offset++;
        this.length = Math.max(this.length, this.offset);
    }

    writeUnsignedShort(value) {
        if (this.readMode) throw new Error("Buffer is in read mode");
        this._ensureCapacity(2);
        this.dataView.setUint16(this.offset, value & 0xFFFF, false); // Big-endian
        this.offset += 2;
        this.length = Math.max(this.length, this.offset);
    }

    // Writes a "smart" unsigned short. If value < 128, writes 1 byte. Else, writes 2 bytes (with MSB set).
    writeUShortSmart(value) {
        if (value < 128) {
            this.writeUnsignedByte(value);
        } else {
            const encodedValue = value + 32768;
             if (encodedValue < 0 || encodedValue > 65535) {
                console.warn("Value out of range for 2-byte UShortSmart encoding: " + value + " (encoded: " + encodedValue + ")");
                // Fallback or error handling might be needed here depending on strictness
                this.writeUnsignedByte(0x80 | ((32768 >>> 8) & 0xFF)); // Write a default smart short like 0
                this.writeUnsignedByte(32768 & 0xFF);
                return;
            }
            this.writeUnsignedByte(((encodedValue >>> 8) & 0xFF) | 0x80); // Set MSB on the first byte
            this.writeUnsignedByte(encodedValue & 0xFF);
        }
    }

    // Writes a "smart" signed short. Handles values in a specific range.
    writeShortSmart(value) {
        if (value >= -64 && value < 64) {
            this.writeUnsignedByte(value + 64);
        } else {
            const unsignedShortValue = value + 49152; // Adjust to unsigned range
            this.writeUnsignedByte(((unsignedShortValue >>> 8) & 0xFF) | 0x80); // Set MSB
            this.writeUnsignedByte(unsignedShortValue & 0xFF);
        }
    }

    writeInt(value) {
        if (this.readMode) throw new Error("Buffer is in read mode");
        this._ensureCapacity(4);
        this.dataView.setInt32(this.offset, value, false); // Big-endian
        this.offset += 4;
        this.length = Math.max(this.length, this.offset);
    }

    readUnsignedByte() {
        if (!this.readMode) throw new Error("Buffer is in write mode");
        if (this.offset >= this.length) throw new Error("Buffer underflow");
        const value = this.dataView.getUint8(this.offset);
        this.offset++;
        return value;
    }

    readUnsignedShort() {
        if (!this.readMode) throw new Error("Buffer is in write mode");
        if (this.offset + 1 >= this.length) throw new Error("Buffer underflow");
        const value = this.dataView.getUint16(this.offset, false); // Big-endian
        this.offset += 2;
        return value;
    }

    readInt() {
        if (!this.readMode) throw new Error("Buffer is in write mode");
        if (this.offset + 3 >= this.length) throw new Error("Buffer underflow");
        const value = this.dataView.getInt32(this.offset, false); // Big-endian
        this.offset += 4;
        return value;
    }

    // Reads a "smart" signed short.
    readShortSmart() {
        const peek = this.dataView.getUint8(this.offset);
        return peek < 128 ? this.readUnsignedByte() - 64 : this.readUnsignedShort() - 49152;
    }

    // Reads a "smart" unsigned short.
    readUShortSmart() {
        const peek = this.dataView.getUint8(this.offset);
        return peek < 128 ? this.readUnsignedByte() : this.readUnsignedShort() - 32768;
    }

    getBytes() {
        if (this.readMode) {
            return new Uint8Array(this.arrayBuffer, 0, this.length);
        }
        // In write mode, return the written portion
        return new Uint8Array(this.arrayBuffer, 0, this.length);
    }
}

// --- SoundEnvelopeJS Class ---
class SoundEnvelopeJS {
    constructor() {
        this.segments = 2;
        this.durations = [0, 65535]; // Durations of each segment
        this.phases = [0, 65535];    // Target phase/value at the end of each segment
        this.start = 0;             // Starting value of the envelope
        this.end = 0;               // Ending value of the envelope
        this.form = 0;              // Shape of the envelope (0: Off, 1: Square, 2: Sine, 3: Triangle, 4: Noise)
        this.ticks = 0;             // Internal counter for envelope progression
        this.phaseIndex = 0;        // Current segment index
        this.step = 0;              // Value change per tick
        this.amplitude = 0;         // Current amplitude/value
        this.max = 0;               // Max ticks for current segment
    }

    decode(buffer) {
        this.form = buffer.readUnsignedByte();
        this.start = buffer.readInt();
        this.end = buffer.readInt();
        this.decodeSegments(buffer);
    }

    decodeSegments(buffer) {
        this.segments = buffer.readUnsignedByte();
        this.durations = new Array(this.segments);
        this.phases = new Array(this.segments);
        for (let i = 0; i < this.segments; ++i) {
            this.durations[i] = buffer.readUnsignedShort();
            this.phases[i] = buffer.readUnsignedShort();
        }
    }

    reset() {
        this.ticks = 0;
        this.phaseIndex = 0;
        this.step = 0;
        this.amplitude = 0;
        this.max = 0;
    }

    doStep(period) {
        if (this.max >= this.ticks) {
            this.amplitude = this.phases[this.phaseIndex++] << 15;
            if (this.phaseIndex >= this.segments) {
                this.phaseIndex = this.segments - 1;
            }
            this.ticks = Math.floor((this.durations[this.phaseIndex] / 65536.0) * period);
            if (this.ticks > this.max) {
                this.step = Math.floor(((this.phases[this.phaseIndex] << 15) - this.amplitude) / (this.ticks - this.max));
            }
        }
        this.amplitude += this.step;
        ++this.max;
        return (this.amplitude - this.step) >> 15;
    }

    encode(buffer) {
        buffer.writeUnsignedByte(this.form);
        buffer.writeInt(this.start);
        buffer.writeInt(this.end);
        this.encodeSegments(buffer);
    }

    encodeSegments(buffer) {
        buffer.writeUnsignedByte(this.segments);
        for (let i = 0; i < this.segments; ++i) {
            buffer.writeUnsignedShort(this.durations[i]);
            buffer.writeUnsignedShort(this.phases[i]);
        }
    }
}

// --- SoundFilterJS Class ---
class SoundFilterJS {
    static minimizedCoefficients = [new Float32Array(8), new Float32Array(8)];
    static coefficients = [new Int32Array(8), new Int32Array(8)];
    static forwardMinimizedCoefficientMultiplier = 0.0;
    static forwardMultiplier = 0;

    constructor() {
        this.pairs = [0, 0]; // Number of pairs for each direction
        this.phases = [
            [new Uint16Array(4), new Uint16Array(4)], // direction 0, phaseType 0 and 1
            [new Uint16Array(4), new Uint16Array(4)]  // direction 1, phaseType 0 and 1
        ];
        this.magnitudes = [
            [new Uint16Array(4), new Uint16Array(4)],
            [new Uint16Array(4), new Uint16Array(4)]
        ];
        this.unity = [0, 0]; // Unity gain values
    }

    adaptMagnitude(direction, i, f) {
        let alpha = this.magnitudes[direction][0][i] + f * (this.magnitudes[direction][1][i] - this.magnitudes[direction][0][i]);
        alpha *= 0.0015258789;
        return 1.0 - Math.pow(10.0, -alpha / 20.0);
    }

    adaptPhase(direction, i, f) {
        let alpha = this.phases[direction][0][i] + f * (this.phases[direction][1][i] - this.phases[direction][0][i]);
        alpha *= 1.2207031E-4;
        return SoundFilterJS.normalize(alpha);
    }

    compute(direction, f) {
        let magnitude;
        if (direction === 0) {
            magnitude = this.unity[0] + (this.unity[1] - this.unity[0]) * f;
            magnitude *= 0.0030517578;
            SoundFilterJS.forwardMinimizedCoefficientMultiplier = Math.pow(0.1, magnitude / 20.0);
            SoundFilterJS.forwardMultiplier = Math.floor(SoundFilterJS.forwardMinimizedCoefficientMultiplier * 65536.0);
        }

        if (this.pairs[direction] === 0) {
            return 0;
        }

        magnitude = this.adaptMagnitude(direction, 0, f);
        SoundFilterJS.minimizedCoefficients[direction][0] = -2.0 * magnitude * Math.cos(this.adaptPhase(direction, 0, f));
        SoundFilterJS.minimizedCoefficients[direction][1] = magnitude * magnitude;

        for (let pair = 1; pair < this.pairs[direction]; ++pair) {
            magnitude = this.adaptMagnitude(direction, pair, f);
            const phase = -2.0 * magnitude * Math.cos(this.adaptPhase(direction, pair, f));
            const coefficient = magnitude * magnitude;
            SoundFilterJS.minimizedCoefficients[direction][pair * 2 + 1] = SoundFilterJS.minimizedCoefficients[direction][pair * 2 - 1] * coefficient;
            SoundFilterJS.minimizedCoefficients[direction][pair * 2] = SoundFilterJS.minimizedCoefficients[direction][pair * 2 - 1] * phase + SoundFilterJS.minimizedCoefficients[direction][pair * 2 - 2] * coefficient;

            for (let pair2 = pair * 2 - 1; pair2 >= 2; --pair2) {
                SoundFilterJS.minimizedCoefficients[direction][pair2] += SoundFilterJS.minimizedCoefficients[direction][pair2 - 1] * phase + SoundFilterJS.minimizedCoefficients[direction][pair2 - 2] * coefficient;
            }
            SoundFilterJS.minimizedCoefficients[direction][1] += SoundFilterJS.minimizedCoefficients[direction][0] * phase + coefficient;
            SoundFilterJS.minimizedCoefficients[direction][0] += phase;
        }

        if (direction === 0) {
            for (let pair = 0; pair < this.pairs[0] * 2; ++pair) {
                SoundFilterJS.minimizedCoefficients[0][pair] *= SoundFilterJS.forwardMinimizedCoefficientMultiplier;
            }
        }

        for (let pair = 0; pair < this.pairs[direction] * 2; ++pair) {
            SoundFilterJS.coefficients[direction][pair] = Math.floor(SoundFilterJS.minimizedCoefficients[direction][pair] * 65536.0);
        }
        return this.pairs[direction] * 2;
    }

    decode(buffer, envelope) {
        const count = buffer.readUnsignedByte();
        this.pairs[0] = count >> 4;
        this.pairs[1] = count & 15;

        if (count !== 0) {
            this.unity[0] = buffer.readUnsignedShort();
            this.unity[1] = buffer.readUnsignedShort();
            const migrated = buffer.readUnsignedByte();

            for (let direction = 0; direction < 2; ++direction) {
                for (let pair = 0; pair < this.pairs[direction]; ++pair) {
                    this.phases[direction][0][pair] = buffer.readUnsignedShort();
                    this.magnitudes[direction][0][pair] = buffer.readUnsignedShort();
                }
            }
            for (let direction = 0; direction < 2; ++direction) {
                for (let pair = 0; pair < this.pairs[direction]; ++pair) {
                    if ((migrated & (1 << (direction * 4 + pair))) !== 0) {
                        this.phases[direction][1][pair] = buffer.readUnsignedShort();
                        this.magnitudes[direction][1][pair] = buffer.readUnsignedShort();
                    } else {
                        this.phases[direction][1][pair] = this.phases[direction][0][pair];
                        this.magnitudes[direction][1][pair] = this.magnitudes[direction][0][pair];
                    }
                }
            }
            if (migrated !== 0 || this.unity[1] !== this.unity[0]) {
                envelope.decodeSegments(buffer);
            }
        } else {
            this.unity[0] = 0;
            this.unity[1] = 0;
        }
    }

    static normalize(alpha) {
        const f = 32.703197 * Math.pow(2.0, alpha);
        return f * Math.PI / 11025.0;
    }

    encode(buffer, envelope) {
        const clampedPairs0 = Math.min(this.pairs[0], 4);
        const clampedPairs1 = Math.min(this.pairs[1], 4);
        const count = (clampedPairs0 << 4) | (clampedPairs1 & 0xF);
        buffer.writeUnsignedByte(count);

        if (count !== 0) {
            buffer.writeUnsignedShort(this.unity[0]);
            buffer.writeUnsignedShort(this.unity[1]);

            const migrated = this._getMigrated(clampedPairs0, clampedPairs1);
            buffer.writeUnsignedByte(migrated);

            for (let direction = 0; direction < 2; ++direction) {
                for (let pair = 0; pair < (direction === 0 ? clampedPairs0 : clampedPairs1); ++pair) {
                    buffer.writeUnsignedShort(this.phases[direction][0][pair]);
                    buffer.writeUnsignedShort(this.magnitudes[direction][0][pair]);
                }
            }
            for (let direction = 0; direction < 2; ++direction) {
                for (let pair = 0; pair < (direction === 0 ? clampedPairs0 : clampedPairs1); ++pair) {
                    if ((migrated & (1 << (direction * 4 + pair))) !== 0) {
                        buffer.writeUnsignedShort(this.phases[direction][1][pair]);
                        buffer.writeUnsignedShort(this.magnitudes[direction][1][pair]);
                    }
                }
            }
            if (migrated !== 0 || this.unity[1] !== this.unity[0]) {
                envelope.encodeSegments(buffer);
            }
        }
    }

    _getMigrated(clampedPairs0, clampedPairs1) {
        let migrated = 0;
        for (let direction = 0; direction < 2; ++direction) {
            for (let pair = 0; pair < (direction === 0 ? clampedPairs0 : clampedPairs1); ++pair) {
                if (this.phases[direction][1][pair] !== this.phases[direction][0][pair] ||
                    this.magnitudes[direction][1][pair] !== this.magnitudes[direction][0][pair]) {
                    migrated |= (1 << (direction * 4 + pair));
                }
            }
        }
        return migrated;
    }
}


// --- SoundToneJS Class ---
class SoundToneJS {
    static toneSamples; // Int32Array
    static toneNoise = new Int32Array(32768);
    static toneSine = new Int32Array(32768);
    static tonePhases = new Int32Array(5);
    static toneDelays = new Int32Array(5);
    static toneVolumeSteps = new Int32Array(5);
    static tonePitchSteps = new Int32Array(5);
    static tonePitchBaseSteps = new Int32Array(5);

    static { // Static initializer block
        for (let i = 0; i < 32768; ++i) {
            SoundToneJS.toneNoise[i] = (Math.floor(Math.random() * 0xFFFFFFFF) & 2) - 1; // Generates -1 or 1
        }
        for (let i = 0; i < 32768; ++i) {
            SoundToneJS.toneSine[i] = Math.floor(Math.sin(i / 5215.1903) * 16384.0);
        }
        // toneSamples will be initialized dynamically based on synthesis needs
    }

    constructor() {
        this.pitch = new SoundEnvelopeJS();
        this.volume = new SoundEnvelopeJS();
        this.pitchModifier = null;
        this.pitchModifierAmplitude = null;
        this.volumeMultiplier = null;
        this.volumeMultiplierAmplitude = null;
        this.release = null;
        this.attack = null;
        this.oscillatorVolume = [0, 0, 0, 0, 0];
        this.oscillatorPitch = [0, 0, 0, 0, 0];
        this.oscillatorDelays = [0, 0, 0, 0, 0];
        this.delayTime = 0;
        this.delayDecay = 100;
        this.filter = new SoundFilterJS();
        this.filterEnvelope = new SoundEnvelopeJS();
        this.duration = 500; // ms
        this.offset = 0;   // ms
    }

    synthesize(steps, tonesDuration) { // tonesDuration is 'tones' in Java, effectively duration in ms
        if (!SoundToneJS.toneSamples || SoundToneJS.toneSamples.length < steps) {
            SoundToneJS.toneSamples = new Int32Array(steps + 22050); // Add some buffer
        }
        // Clear sample buffer
        for (let i = 0; i < steps; i++) SoundToneJS.toneSamples[i] = 0;


        if (tonesDuration >= 10) { // Java: if (tones >= 10)
            const durationFactor = steps / (tonesDuration + 0.0);

            this.pitch.reset();
            this.volume.reset();
            let pitchModulationStep = 0;
            let pitchModulationBaseStep = 0;
            let pitchModulationPhase = 0;
            if (this.pitchModifier) {
                this.pitchModifier.reset();
                this.pitchModifierAmplitude.reset();
                pitchModulationStep = Math.floor((this.pitchModifier.end - this.pitchModifier.start) * 32.768 / durationFactor);
                pitchModulationBaseStep = Math.floor(this.pitchModifier.start * 32.768 / durationFactor);
            }

            let volumeModulationStep = 0;
            let volumeModulationBaseStep = 0;
            let volumeModulationPhase = 0;
            if (this.volumeMultiplier) {
                this.volumeMultiplier.reset();
                this.volumeMultiplierAmplitude.reset();
                volumeModulationStep = Math.floor((this.volumeMultiplier.end - this.volumeMultiplier.start) * 32.768 / durationFactor);
                volumeModulationBaseStep = Math.floor(this.volumeMultiplier.start * 32.768 / durationFactor);
            }

            for (let i = 0; i < 5; ++i) {
                if (this.oscillatorVolume[i] !== 0) {
                    SoundToneJS.tonePhases[i] = 0;
                    SoundToneJS.toneDelays[i] = Math.floor(this.oscillatorDelays[i] * durationFactor);
                    SoundToneJS.toneVolumeSteps[i] = (this.oscillatorVolume[i] << 14) / 100;
                    SoundToneJS.tonePitchSteps[i] = Math.floor((this.pitch.end - this.pitch.start) * 32.768 * Math.pow(1.0057929410678534, this.oscillatorPitch[i]) / durationFactor);
                    SoundToneJS.tonePitchBaseSteps[i] = Math.floor(this.pitch.start * 32.768 / durationFactor);
                }
            }

            let pitchChange, volumeChange, volMulChange, volMulAmpChange;
            for (let step = 0; step < steps; ++step) {
                pitchChange = this.pitch.doStep(steps);
                volumeChange = this.volume.doStep(steps);

                if (this.pitchModifier) {
                    volMulChange = this.pitchModifier.doStep(steps);
                    volMulAmpChange = this.pitchModifierAmplitude.doStep(steps);
                    pitchChange += this.evaluateWave(pitchModulationPhase, volMulAmpChange, this.pitchModifier.form) >> 1;
                    pitchModulationPhase += pitchModulationBaseStep + (volMulChange * pitchModulationStep >> 16);
                }

                if (this.volumeMultiplier) {
                    volMulChange = this.volumeMultiplier.doStep(steps);
                    volMulAmpChange = this.volumeMultiplierAmplitude.doStep(steps);
                    volumeChange = volumeChange * ((this.evaluateWave(volumeModulationPhase, volMulAmpChange, this.volumeMultiplier.form) >> 1) + 32768) >> 15;
                    volumeModulationPhase += volumeModulationBaseStep + (volMulChange * volumeModulationStep >> 16);
                }

                for (let osc = 0; osc < 5; ++osc) {
                    if (this.oscillatorVolume[osc] !== 0) {
                        const targetIdx = SoundToneJS.toneDelays[osc] + step;
                        if (targetIdx < steps) {
                            SoundToneJS.toneSamples[targetIdx] += this.evaluateWave(SoundToneJS.tonePhases[osc], volumeChange * SoundToneJS.toneVolumeSteps[osc] >> 15, this.pitch.form);
                            SoundToneJS.tonePhases[osc] += (pitchChange * SoundToneJS.tonePitchSteps[osc] >> 16) + SoundToneJS.tonePitchBaseSteps[osc];
                        }
                    }
                }
            }

            if (this.release) {
                this.release.reset();
                this.attack.reset();
                let tickCounter = 0;
                let muted = true;
                for (let i = 0; i < steps; ++i) {
                    const releaseVal = this.release.doStep(steps);
                    const attackVal = this.attack.doStep(steps);
                    let threshold;
                    if (muted) {
                        threshold = (releaseVal * (this.release.end - this.release.start) >> 8) + this.release.start;
                    } else {
                        threshold = (attackVal * (this.release.end - this.release.start) >> 8) + this.release.start; // Java uses release.end - release.start for both
                    }
                    tickCounter += 256;
                    if (tickCounter >= threshold) {
                        tickCounter = 0;
                        muted = !muted;
                    }
                    if (muted) {
                        SoundToneJS.toneSamples[i] = 0;
                    }
                }
            }

            if (this.delayTime > 0 && this.delayDecay > 0) {
                const delaySamples = Math.floor(this.delayTime * durationFactor);
                for (let i = delaySamples; i < steps; ++i) {
                    SoundToneJS.toneSamples[i] += SoundToneJS.toneSamples[i - delaySamples] * this.delayDecay / 100;
                }
            }

            if (this.filter.pairs[0] > 0 || this.filter.pairs[1] > 0) {
                this.filterEnvelope.reset();
                let filterStepVal = this.filterEnvelope.doStep(steps + 1);
                let inputCoeffCount = this.filter.compute(0, filterStepVal / 65536.0);
                let outputCoeffCount = this.filter.compute(1, filterStepVal / 65536.0);

                if (steps >= inputCoeffCount + outputCoeffCount) {
                    let currentSampleIdx = 0;
                    let limit1 = Math.min(outputCoeffCount, steps - inputCoeffCount); // Java: volumeMultiplierAmplitudeChange

                    while(currentSampleIdx < limit1) {
                        let filteredSample = (SoundToneJS.toneSamples[currentSampleIdx + inputCoeffCount] * SoundFilterJS.forwardMultiplier) >> 16;
                        for(let k=0; k < inputCoeffCount; ++k) {
                            filteredSample += (SoundToneJS.toneSamples[currentSampleIdx + inputCoeffCount - 1 - k] * SoundFilterJS.coefficients[0][k]) >> 16;
                        }
                        for(let k=0; k < currentSampleIdx; ++k) { // Java: var17 < volumeMultiplierChange (currentSampleIdx)
                            filteredSample -= (SoundToneJS.toneSamples[currentSampleIdx - 1 - k] * SoundFilterJS.coefficients[1][k]) >> 16;
                        }
                        SoundToneJS.toneSamples[currentSampleIdx] = filteredSample;
                        filterStepVal = this.filterEnvelope.doStep(steps + 1);
                        ++currentSampleIdx;
                    }

                    let limit2Base = 128; // Java: volumeMultiplierAmplitudeChange
                    while(true) {
                        let limit2 = Math.min(limit2Base, steps - inputCoeffCount);
                        if (currentSampleIdx >= limit2) limit2 = steps - inputCoeffCount; // Ensure we process remaining samples

                        while(currentSampleIdx < limit2) {
                            let filteredSample = (SoundToneJS.toneSamples[currentSampleIdx + inputCoeffCount] * SoundFilterJS.forwardMultiplier) >> 16;
                            for(let k=0; k < inputCoeffCount; ++k) {
                                filteredSample += (SoundToneJS.toneSamples[currentSampleIdx + inputCoeffCount - 1 - k] * SoundFilterJS.coefficients[0][k]) >> 16;
                            }
                            for(let k=0; k < outputCoeffCount; ++k) { // Java: var18 < volumeChange (outputCoeffCount)
                                filteredSample -= (SoundToneJS.toneSamples[currentSampleIdx - 1 - k] * SoundFilterJS.coefficients[1][k]) >> 16;
                            }
                            SoundToneJS.toneSamples[currentSampleIdx] = filteredSample;
                            filterStepVal = this.filterEnvelope.doStep(steps + 1);
                            ++currentSampleIdx;
                        }

                        if (currentSampleIdx >= steps - inputCoeffCount) {
                            while(currentSampleIdx < steps) {
                                let filteredSample = 0;
                                for(let k = currentSampleIdx + inputCoeffCount - steps; k < inputCoeffCount; ++k) {
                                    filteredSample += (SoundToneJS.toneSamples[currentSampleIdx + inputCoeffCount - 1 - k] * SoundFilterJS.coefficients[0][k]) >> 16;
                                }
                                for(let k=0; k < outputCoeffCount; ++k) {
                                    filteredSample -= (SoundToneJS.toneSamples[currentSampleIdx - 1 - k] * SoundFilterJS.coefficients[1][k]) >> 16;
                                }
                                SoundToneJS.toneSamples[currentSampleIdx] = filteredSample;
                                this.filterEnvelope.doStep(steps + 1); // filterStepVal not used here
                                ++currentSampleIdx;
                            }
                            break;
                        }
                        inputCoeffCount = this.filter.compute(0, filterStepVal / 65536.0);
                        outputCoeffCount = this.filter.compute(1, filterStepVal / 65536.0);
                        limit2Base += 128;
                    }
                }
            }


            for (let i = 0; i < steps; ++i) {
                if (SoundToneJS.toneSamples[i] < -32768) SoundToneJS.toneSamples[i] = -32768;
                if (SoundToneJS.toneSamples[i] > 32767) SoundToneJS.toneSamples[i] = 32767;
            }
        }
        // Return a slice of the samples up to 'steps' length
        return SoundToneJS.toneSamples.slice(0, steps);
    }


    evaluateWave(phase, amplitude, form) {
        if (form === 1) { // Square
            return (phase & 32767) < 16384 ? amplitude : -amplitude;
        } else if (form === 2) { // Sine
            return SoundToneJS.toneSine[phase & 32767] * amplitude >> 14;
        } else if (form === 3) { // Triangle
            return (amplitude * (phase & 32767) >> 14) - amplitude;
        } else if (form === 4) { // Noise
            return amplitude * SoundToneJS.toneNoise[Math.floor(phase / 2607) & 32767];
        }
        return 0; // Off
    }

    decode(buffer) {
        this.pitch = new SoundEnvelopeJS(); this.pitch.decode(buffer);
        this.volume = new SoundEnvelopeJS(); this.volume.decode(buffer);

        let flag = buffer.readUnsignedByte();
        if (flag !== 0) {
            buffer.offset--;
            this.pitchModifier = new SoundEnvelopeJS(); this.pitchModifier.decode(buffer);
            this.pitchModifierAmplitude = new SoundEnvelopeJS(); this.pitchModifierAmplitude.decode(buffer);
        } else { this.pitchModifier = null; this.pitchModifierAmplitude = null; }

        flag = buffer.readUnsignedByte();
        if (flag !== 0) {
            buffer.offset--;
            this.volumeMultiplier = new SoundEnvelopeJS(); this.volumeMultiplier.decode(buffer);
            this.volumeMultiplierAmplitude = new SoundEnvelopeJS(); this.volumeMultiplierAmplitude.decode(buffer);
        } else { this.volumeMultiplier = null; this.volumeMultiplierAmplitude = null; }

        flag = buffer.readUnsignedByte();
        if (flag !== 0) {
            buffer.offset--;
            this.release = new SoundEnvelopeJS(); this.release.decode(buffer);
            this.attack = new SoundEnvelopeJS(); this.attack.decode(buffer);
        } else { this.release = null; this.attack = null; }

        for (let i = 0; i < 5; ++i) {
            const vol = buffer.readUShortSmart();
            if (vol === 0) {
                 for (let j = i; j < 5; j++) {
                    this.oscillatorVolume[j] = 0;
                    this.oscillatorPitch[j] = 0;
                    this.oscillatorDelays[j] = 0;
                }
                break;
            }
            this.oscillatorVolume[i] = vol;
            this.oscillatorPitch[i] = buffer.readShortSmart();
            this.oscillatorDelays[i] = buffer.readUShortSmart();
        }

        this.delayTime = buffer.readUShortSmart();
        this.delayDecay = buffer.readUShortSmart();
        this.duration = buffer.readUnsignedShort();
        this.offset = buffer.readUnsignedShort();

        this.filter = new SoundFilterJS();
        this.filterEnvelope = new SoundEnvelopeJS();
        this.filter.decode(buffer, this.filterEnvelope);
    }

    encode(buffer) {
        this.pitch.encode(buffer);
        this.volume.encode(buffer);

        if (this.pitchModifier && this.pitchModifier.form !== 0) {
            this.pitchModifier.encode(buffer);
            if (this.pitchModifierAmplitude) this.pitchModifierAmplitude.encode(buffer); else new SoundEnvelopeJS().encode(buffer); // Ensure something is written
        } else { buffer.writeUnsignedByte(0); }

        if (this.volumeMultiplier && this.volumeMultiplier.form !== 0) {
            this.volumeMultiplier.encode(buffer);
             if (this.volumeMultiplierAmplitude) this.volumeMultiplierAmplitude.encode(buffer); else new SoundEnvelopeJS().encode(buffer);
        } else { buffer.writeUnsignedByte(0); }

        if (this.release && this.release.form !== 0) {
            this.release.encode(buffer);
            if (this.attack) this.attack.encode(buffer); else new SoundEnvelopeJS().encode(buffer);
        } else { buffer.writeUnsignedByte(0); }

        let terminatorWritten = false;
        for (let i = 0; i < 5; ++i) {
            if (this.oscillatorVolume[i] !== 0 || this.oscillatorPitch[i] !== 0 || this.oscillatorDelays[i] !== 0) {
                buffer.writeUShortSmart(this.oscillatorVolume[i]);
                buffer.writeShortSmart(this.oscillatorPitch[i]);
                buffer.writeUShortSmart(this.oscillatorDelays[i]);
            } else {
                buffer.writeUShortSmart(0); // Terminator
                terminatorWritten = true;
                break;
            }
        }
         if (!terminatorWritten) { // If all 5 oscillators were used, still need a terminator
            buffer.writeUShortSmart(0);
        }


        buffer.writeUShortSmart(this.delayTime);
        buffer.writeUShortSmart(this.delayDecay);
        buffer.writeUnsignedShort(this.duration);
        buffer.writeUnsignedShort(this.offset);
        this.filter.encode(buffer, this.filterEnvelope);
    }
}

// --- SoundEffectJS Class ---
class SoundEffectJS {
    constructor(buffer) {
        this.soundTones = new Array(10).fill(null);
        this.start = 0;
        this.end = 0;

        if (buffer instanceof BufferJS) {
            for (let i = 0; i < 10; ++i) {
                const toneExistsFlag = buffer.readUnsignedByte();
                if (toneExistsFlag !== 0) {
                    buffer.offset--; // Rewind
                    this.soundTones[i] = new SoundToneJS();
                    this.soundTones[i].decode(buffer);
                }
            }
            this.start = buffer.readUnsignedShort();
            this.end = buffer.readUnsignedShort();
        }
    }

    mix() { // Returns 8-bit signed byte array
        let maxTotalDurationMs = 0;
        for (let i = 0; i < 10; ++i) {
            if (this.soundTones[i] && (this.soundTones[i].duration + this.soundTones[i].offset > maxTotalDurationMs)) {
                maxTotalDurationMs = this.soundTones[i].duration + this.soundTones[i].offset;
            }
        }

        if (maxTotalDurationMs === 0) return new Int8Array(0);

        const totalSamples = Math.floor(maxTotalDurationMs * 22050 / 1000);
        const mixedSamples8Bit = new Int8Array(totalSamples); // For 8-bit output

        for (let i = 0; i < 10; ++i) {
            const tone = this.soundTones[i];
            if (tone) {
                const toneDurationSamples = Math.floor(tone.duration * 22050 / 1000);
                const toneOffsetSamples = Math.floor(tone.offset * 22050 / 1000);
                const synthesized16Bit = tone.synthesize(toneDurationSamples, tone.duration); // Pass duration in ms

                for (let j = 0; j < toneDurationSamples; ++j) {
                    if (toneOffsetSamples + j < totalSamples) {
                        let sample16 = synthesized16Bit[j];
                        // Convert 16-bit to 8-bit and mix
                        let sample8 = sample16 >> 8;
                        let currentMix = mixedSamples8Bit[toneOffsetSamples + j] + sample8;

                        // Clamp to 8-bit signed range (-128 to 127)
                        if (currentMix < -128) currentMix = -128;
                        else if (currentMix > 127) currentMix = 127;
                        mixedSamples8Bit[toneOffsetSamples + j] = currentMix;
                    }
                }
            }
        }
        return mixedSamples8Bit;
    }

    encode() {
        const buffer = new BufferJS(); // Write mode
        for (let i = 0; i < 10; ++i) {
            if (this.soundTones[i]) {
                this.soundTones[i].encode(buffer);
            } else {
                buffer.writeUnsignedByte(0);
            }
        }
        buffer.writeUnsignedShort(this.start);
        buffer.writeUnsignedShort(this.end);
        return buffer.getBytes();
    }
}


// --- SoundEffectEditorJS (Main Application Logic) ---
class SoundEffectEditorJS {
    constructor() {
        this.audioContext = new (window.AudioContext || window.webkitAudioContext)();
        this.currentSoundEffect = new SoundEffectJS(); // Start with an empty effect
        this.currentSoundEffect.soundTones[0] = new SoundToneJS(); // Add one default tone
        this.selectedToneIndex = 0;
        this.selectedOscillatorIndex = 0;
        this.currentAudioBufferSource = null;
        this.outputBitDepth = 16; // Default

        this.ui = {
            importFile: document.getElementById('importFile'),
            saveJagFX: document.getElementById('saveJagFX'),
            exportWav: document.getElementById('exportWav'),
            outputBitDepthSelect: document.getElementById('outputBitDepth'),
            playFullEffect: document.getElementById('playFullEffect'),
            playSelectedTone: document.getElementById('playSelectedTone'),
            stopSound: document.getElementById('stopSound'),
            addTone: document.getElementById('addTone'),
            removeTone: document.getElementById('removeTone'),
            toneListDiv: document.getElementById('toneList'),
            selectedToneControlsDiv: document.getElementById('selectedToneControls'),
            currentToneIndexDisplay: document.getElementById('currentToneIndexDisplay'),
            durationInput: document.getElementById('duration'),
            offsetInput: document.getElementById('offset'),
            delayTimeInput: document.getElementById('delayTime'),
            delayDecayInput: document.getElementById('delayDecay'),
            oscillatorSelect: document.getElementById('oscillatorSelect'),
            currentOscillatorIndexDisplay: document.getElementById('currentOscillatorIndexDisplay'),
            oscVolumeInput: document.getElementById('oscVolume'),
            oscPitchInput: document.getElementById('oscPitch'),
            oscDelayInput: document.getElementById('oscDelay'),
            // Envelope UI
            pitchEnvForm: document.getElementById('pitchEnvForm'),
            pitchEnvStart: document.getElementById('pitchEnvStart'),
            pitchEnvEnd: document.getElementById('pitchEnvEnd'),
            volumeEnvForm: document.getElementById('volumeEnvForm'),
            volumeEnvStart: document.getElementById('volumeEnvStart'),
            volumeEnvEnd: document.getElementById('volumeEnvEnd'),
            // Optional Envelopes
            pitchModifierEnabledCheckbox: document.getElementById('pitchModifierEnabled'),
            pitchModifierControlsDiv: document.getElementById('pitchModifierControls'),
            pitchModEnvForm: document.getElementById('pitchModEnvForm'),
            pitchModEnvStart: document.getElementById('pitchModEnvStart'),
            pitchModEnvEnd: document.getElementById('pitchModEnvEnd'),
            pitchModAmpEnvForm: document.getElementById('pitchModAmpEnvForm'),
            pitchModAmpEnvStart: document.getElementById('pitchModAmpEnvStart'),
            pitchModAmpEnvEnd: document.getElementById('pitchModAmpEnvEnd'),
            // Filter UI
            filterEnabledCheckbox: document.getElementById('filterEnabled'),
            filterControlsDiv: document.getElementById('filterControls'),
            filterUnity0Input: document.getElementById('filterUnity0'),
            filterUnity1Input: document.getElementById('filterUnity1'),
            filterPairs0Input: document.getElementById('filterPairs0'),
            filterPairs1Input: document.getElementById('filterPairs1'),
            filterCoefficientsPanel: document.getElementById('filterCoefficientsPanel'),


            waveformCanvas: document.getElementById('waveformCanvas'),
            pitchEnvelopeCanvas: document.getElementById('pitchEnvelopeCanvas'),
            volumeEnvelopeCanvas: document.getElementById('volumeEnvelopeCanvas'),
            filterCanvas: document.getElementById('filterCanvas'),
        };

        this.waveformCtx = this.ui.waveformCanvas.getContext('2d');
        this.pitchEnvCtx = this.ui.pitchEnvelopeCanvas.getContext('2d');
        this.volumeEnvCtx = this.ui.volumeEnvelopeCanvas.getContext('2d');
        this.filterCtx = this.ui.filterCanvas.getContext('2d');

        this.debounceTimer = null;

        this.initEventListeners();
        this.renderToneList();
        this.selectTone(0);
    }

    _clearDebounce() {
        if (this.debounceTimer) clearTimeout(this.debounceTimer);
    }

    _triggerDebounceUpdate() {
        this._clearDebounce();
        this.debounceTimer = setTimeout(() => {
            this.updateFromUI(); // Apply changes to model
            this.drawVisualizations(); // Redraw graphs
        }, 250); // 250ms debounce
    }


    initEventListeners() {
        this.ui.importFile.addEventListener('change', (e) => this.importJagFX(e));
        this.ui.saveJagFX.addEventListener('click', () => this.saveJagFX());
        this.ui.exportWav.addEventListener('click', () => this.exportToWav());
        this.ui.outputBitDepthSelect.addEventListener('change', (e) => {
            this.outputBitDepth = parseInt(e.target.value);
            this.drawVisualizations(); // Redraw waveform with new bit depth scaling
        });

        this.ui.playFullEffect.addEventListener('click', () => this.playFullEffect());
        this.ui.playSelectedTone.addEventListener('click', () => this.playSelectedTone());
        this.ui.stopSound.addEventListener('click', () => this.stopSound());
        this.ui.addTone.addEventListener('click', () => this.addTone());
        this.ui.removeTone.addEventListener('click', () => this.removeTone());

        // Tone parameter inputs - use 'input' for live updates, or 'change'
        const inputsToListen = [
            this.ui.durationInput, this.ui.offsetInput, this.ui.delayTimeInput, this.ui.delayDecayInput,
            this.ui.oscVolumeInput, this.ui.oscPitchInput, this.ui.oscDelayInput,
            this.ui.pitchEnvForm, this.ui.pitchEnvStart, this.ui.pitchEnvEnd,
            this.ui.volumeEnvForm, this.ui.volumeEnvStart, this.ui.volumeEnvEnd,
            this.ui.pitchModEnvForm, this.ui.pitchModEnvStart, this.ui.pitchModEnvEnd,
            this.ui.pitchModAmpEnvForm, this.ui.pitchModAmpEnvStart, this.ui.pitchModAmpEnvEnd,
            this.ui.filterUnity0Input, this.ui.filterUnity1Input,
            this.ui.filterPairs0Input, this.ui.filterPairs1Input
        ];
        inputsToListen.forEach(input => {
            if (input) input.addEventListener('input', () => this._triggerDebounceUpdate());
        });

        this.ui.oscillatorSelect.addEventListener('change', (e) => {
            this.selectedOscillatorIndex = parseInt(e.target.value);
            this.ui.currentOscillatorIndexDisplay.textContent = this.selectedOscillatorIndex;
            this.updateUIFromModel(); // Refresh oscillator inputs
        });

        this.ui.pitchModifierEnabledCheckbox.addEventListener('change', (e) => {
            this.ui.pitchModifierControlsDiv.classList.toggle('hidden', !e.target.checked);
            this._triggerDebounceUpdate();
        });
        this.ui.filterEnabledCheckbox.addEventListener('change', (e) => {
            this.ui.filterControlsDiv.classList.toggle('hidden', !e.target.checked);
            this._triggerDebounceUpdate();
        });

        // Listener for filter pairs change to update coefficient UI
        this.ui.filterPairs0Input.addEventListener('input', () => this.renderFilterCoefficientsUI());
        this.ui.filterPairs1Input.addEventListener('input', () => this.renderFilterCoefficientsUI());

    }

    getSelectedTone() {
        if (this.selectedToneIndex >= 0 && this.selectedToneIndex < this.currentSoundEffect.soundTones.length) {
            return this.currentSoundEffect.soundTones[this.selectedToneIndex];
        }
        return null;
    }

    updateFromUI() {
        const tone = this.getSelectedTone();
        if (!tone) return;

        tone.duration = parseInt(this.ui.durationInput.value) || 0;
        tone.offset = parseInt(this.ui.offsetInput.value) || 0;
        tone.delayTime = parseInt(this.ui.delayTimeInput.value) || 0;
        tone.delayDecay = parseInt(this.ui.delayDecayInput.value) || 0;

        tone.oscillatorVolume[this.selectedOscillatorIndex] = parseInt(this.ui.oscVolumeInput.value) || 0;
        tone.oscillatorPitch[this.selectedOscillatorIndex] = parseInt(this.ui.oscPitchInput.value) || 0;
        tone.oscillatorDelays[this.selectedOscillatorIndex] = parseInt(this.ui.oscDelayInput.value) || 0;

        // Pitch Envelope
        tone.pitch.form = parseInt(this.ui.pitchEnvForm.value);
        tone.pitch.start = parseInt(this.ui.pitchEnvStart.value) || 0;
        tone.pitch.end = parseInt(this.ui.pitchEnvEnd.value) || 0;

        // Volume Envelope
        tone.volume.form = parseInt(this.ui.volumeEnvForm.value);
        tone.volume.start = parseInt(this.ui.volumeEnvStart.value) || 0;
        tone.volume.end = parseInt(this.ui.volumeEnvEnd.value) || 0;

        // Optional: Pitch Modifier
        if (this.ui.pitchModifierEnabledCheckbox.checked) {
            if (!tone.pitchModifier) {
                tone.pitchModifier = new SoundEnvelopeJS();
                tone.pitchModifierAmplitude = new SoundEnvelopeJS();
            }
            tone.pitchModifier.form = parseInt(this.ui.pitchModEnvForm.value);
            tone.pitchModifier.start = parseInt(this.ui.pitchModEnvStart.value) || 0;
            tone.pitchModifier.end = parseInt(this.ui.pitchModEnvEnd.value) || 0;
            tone.pitchModifierAmplitude.form = parseInt(this.ui.pitchModAmpEnvForm.value);
            tone.pitchModifierAmplitude.start = parseInt(this.ui.pitchModAmpEnvStart.value) || 0;
            tone.pitchModifierAmplitude.end = parseInt(this.ui.pitchModAmpEnvEnd.value) || 0;
        } else {
            tone.pitchModifier = null;
            tone.pitchModifierAmplitude = null;
        }

        // Filter
        if (this.ui.filterEnabledCheckbox.checked) {
            if (!tone.filter) {
                tone.filter = new SoundFilterJS();
                tone.filterEnvelope = new SoundEnvelopeJS(); // Ensure filter envelope exists
            }
            tone.filter.unity[0] = parseInt(this.ui.filterUnity0Input.value) || 0;
            tone.filter.unity[1] = parseInt(this.ui.filterUnity1Input.value) || 0;
            tone.filter.pairs[0] = parseInt(this.ui.filterPairs0Input.value) || 0;
            tone.filter.pairs[1] = parseInt(this.ui.filterPairs1Input.value) || 0;

            // Update filter coefficients from dynamic inputs
            for (let dir = 0; dir < 2; dir++) {
                const numPairs = tone.filter.pairs[dir];
                for (let pair = 0; pair < numPairs; pair++) {
                    const phase0Input = document.getElementById(`filter_dir${dir}_pair${pair}_phase0`);
                    const mag0Input = document.getElementById(`filter_dir${dir}_pair${pair}_mag0`);
                    const phase1Input = document.getElementById(`filter_dir${dir}_pair${pair}_phase1`);
                    const mag1Input = document.getElementById(`filter_dir${dir}_pair${pair}_mag1`);

                    if (phase0Input) tone.filter.phases[dir][0][pair] = parseInt(phase0Input.value) || 0;
                    if (mag0Input) tone.filter.magnitudes[dir][0][pair] = parseInt(mag0Input.value) || 0;
                    if (phase1Input) tone.filter.phases[dir][1][pair] = parseInt(phase1Input.value) || 0;
                    if (mag1Input) tone.filter.magnitudes[dir][1][pair] = parseInt(mag1Input.value) || 0;
                }
            }

        } else {
            tone.filter = null; // Or reset its pairs to 0 if it should always exist
            tone.filterEnvelope = null; // Or reset its form to 0
        }

        // console.log("Updated tone:", tone);
    }

    updateUIFromModel() {
        const tone = this.getSelectedTone();
        this.ui.selectedToneControlsDiv.classList.toggle('hidden', !tone);
        if (!tone) return;

        this.ui.currentToneIndexDisplay.textContent = `Tone ${this.selectedToneIndex + 1}`;
        this.ui.durationInput.value = tone.duration;
        this.ui.offsetInput.value = tone.offset;
        this.ui.delayTimeInput.value = tone.delayTime;
        this.ui.delayDecayInput.value = tone.delayDecay;

        this.ui.oscillatorSelect.value = this.selectedOscillatorIndex;
        this.ui.currentOscillatorIndexDisplay.textContent = this.selectedOscillatorIndex;
        this.ui.oscVolumeInput.value = tone.oscillatorVolume[this.selectedOscillatorIndex];
        this.ui.oscPitchInput.value = tone.oscillatorPitch[this.selectedOscillatorIndex];
        this.ui.oscDelayInput.value = tone.oscillatorDelays[this.selectedOscillatorIndex];

        // Pitch Envelope
        this.ui.pitchEnvForm.value = tone.pitch.form;
        this.ui.pitchEnvStart.value = tone.pitch.start;
        this.ui.pitchEnvEnd.value = tone.pitch.end;

        // Volume Envelope
        this.ui.volumeEnvForm.value = tone.volume.form;
        this.ui.volumeEnvStart.value = tone.volume.start;
        this.ui.volumeEnvEnd.value = tone.volume.end;

        // Optional: Pitch Modifier
        const pmEnabled = !!tone.pitchModifier;
        this.ui.pitchModifierEnabledCheckbox.checked = pmEnabled;
        this.ui.pitchModifierControlsDiv.classList.toggle('hidden', !pmEnabled);
        if (pmEnabled) {
            this.ui.pitchModEnvForm.value = tone.pitchModifier.form;
            this.ui.pitchModEnvStart.value = tone.pitchModifier.start;
            this.ui.pitchModEnvEnd.value = tone.pitchModifier.end;
            this.ui.pitchModAmpEnvForm.value = tone.pitchModifierAmplitude.form;
            this.ui.pitchModAmpEnvStart.value = tone.pitchModifierAmplitude.start;
            this.ui.pitchModAmpEnvEnd.value = tone.pitchModifierAmplitude.end;
        }

        // Filter
        const filterEnabled = !!tone.filter && (tone.filter.pairs[0] > 0 || tone.filter.pairs[1] > 0);
        this.ui.filterEnabledCheckbox.checked = filterEnabled;
        this.ui.filterControlsDiv.classList.toggle('hidden', !filterEnabled);
        if (tone.filter) {
            this.ui.filterUnity0Input.value = tone.filter.unity[0];
            this.ui.filterUnity1Input.value = tone.filter.unity[1];
            this.ui.filterPairs0Input.value = tone.filter.pairs[0];
            this.ui.filterPairs1Input.value = tone.filter.pairs[1];
            this.renderFilterCoefficientsUI();
        } else {
             this.ui.filterUnity0Input.value = 0;
            this.ui.filterUnity1Input.value = 0;
            this.ui.filterPairs0Input.value = 0;
            this.ui.filterPairs1Input.value = 0;
            this.renderFilterCoefficientsUI(); // Clear it
        }


        this.drawVisualizations();
    }

    renderFilterCoefficientsUI() {
        this.ui.filterCoefficientsPanel.innerHTML = ''; // Clear previous
        const tone = this.getSelectedTone();
        if (!tone || !tone.filter || !this.ui.filterEnabledCheckbox.checked) return;

        for (let dir = 0; dir < 2; dir++) {
            const numPairs = parseInt(dir === 0 ? this.ui.filterPairs0Input.value : this.ui.filterPairs1Input.value) || 0;
            if (numPairs > 0) {
                 const dirLabel = document.createElement('h5');
                 dirLabel.textContent = `Direction ${dir} Coefficients (${numPairs} pairs):`;
                 this.ui.filterCoefficientsPanel.appendChild(dirLabel);
            }
            for (let pair = 0; pair < numPairs; pair++) {
                const pairDiv = document.createElement('div');
                pairDiv.className = 'filter-pair-control';
                pairDiv.innerHTML = `
                    <label>P${pair}:</label>
                    <label>Ph0:</label><input type="number" id="filter_dir${dir}_pair${pair}_phase0" value="${tone.filter.phases[dir][0][pair] || 0}">
                    <label>Mg0:</label><input type="number" id="filter_dir${dir}_pair${pair}_mag0" value="${tone.filter.magnitudes[dir][0][pair] || 0}">
                    <label>Ph1:</label><input type="number" id="filter_dir${dir}_pair${pair}_phase1" value="${tone.filter.phases[dir][1][pair] || 0}">
                    <label>Mg1:</label><input type="number" id="filter_dir${dir}_pair${pair}_mag1" value="${tone.filter.magnitudes[dir][1][pair] || 0}">
                `;
                this.ui.filterCoefficientsPanel.appendChild(pairDiv);
                // Add event listeners to these new inputs
                ['phase0', 'mag0', 'phase1', 'mag1'].forEach(type => {
                    document.getElementById(`filter_dir${dir}_pair${pair}_${type}`).addEventListener('input', () => this._triggerDebounceUpdate());
                });
            }
        }
    }


    renderToneList() {
        this.ui.toneListDiv.innerHTML = '';
        this.currentSoundEffect.soundTones.forEach((tone, index) => {
            if (tone) { // Only render existing tones
                const item = document.createElement('div');
                item.textContent = `Tone ${index + 1}`;
                item.className = 'tone-list-item';
                if (index === this.selectedToneIndex) {
                    item.classList.add('selected');
                }
                item.addEventListener('click', () => this.selectTone(index));
                this.ui.toneListDiv.appendChild(item);
            }
        });
         // Ensure at least one tone is selected if list is not empty
        if (this.currentSoundEffect.soundTones.filter(t => t !== null).length > 0 && this.selectedToneIndex === -1) {
            this.selectTone(0);
        } else if (this.currentSoundEffect.soundTones.filter(t => t !== null).length === 0) {
            this.selectTone(-1); // No tones, deselect
        }
    }

    selectTone(index) {
        this.selectedToneIndex = index;
         if (index === -1 || !this.currentSoundEffect.soundTones[index]) {
            this.ui.selectedToneControlsDiv.classList.add('hidden');
        } else {
            this.ui.selectedToneControlsDiv.classList.remove('hidden');
            this.updateUIFromModel();
        }
        this.renderToneList(); // Re-render to update selection highlight
    }

    addTone() {
        if (this.currentSoundEffect.soundTones.filter(t => t !== null).length < 10) {
            const newTone = new SoundToneJS();
            // Find first null spot or append
            let added = false;
            for(let i=0; i < 10; i++) {
                if(this.currentSoundEffect.soundTones[i] === null) {
                    this.currentSoundEffect.soundTones[i] = newTone;
                    this.selectTone(i);
                    added = true;
                    break;
                }
            }
            // This logic might be flawed if we strictly adhere to 10 slots and some are null in between
            // For simplicity, let's assume we just manage a list that can grow up to 10.
            // The Java code always has 10 slots, some might be null.
            // A more robust way for JS might be to filter out nulls then add if length < 10
            // For now, this matches the fixed-size array concept somewhat.
            this.renderToneList();
        } else {
            alert("Maximum of 10 tones reached.");
        }
    }

    removeTone() {
        if (this.selectedToneIndex !== -1 && this.currentSoundEffect.soundTones[this.selectedToneIndex]) {
            if (this.currentSoundEffect.soundTones.filter(t => t !== null).length > 1) {
                this.currentSoundEffect.soundTones[this.selectedToneIndex] = null; // Mark as null
                this.renderToneList();
                // Select the next available tone or the previous one
                let nextSelection = -1;
                for(let i=0; i < 10; i++) { if(this.currentSoundEffect.soundTones[i]) { nextSelection = i; break; } }
                this.selectTone(nextSelection);

            } else {
                alert("Cannot remove the last tone.");
            }
        }
    }

    importJagFX(event) {
        const file = event.target.files[0];
        if (!file) return;
        const reader = new FileReader();
        reader.onload = (e) => {
            try {
                const buffer = new BufferJS(new Uint8Array(e.target.result));
                this.currentSoundEffect = new SoundEffectJS(buffer);
                // Ensure there's at least one tone for the UI, even if file is empty/corrupt
                if (this.currentSoundEffect.soundTones.filter(t => t !== null).length === 0) {
                   this.currentSoundEffect.soundTones[0] = new SoundToneJS();
                }
                this.selectedToneIndex = -1; // Reset selection
                this.renderToneList(); // This will auto-select the first valid tone
                alert('JagFX file imported successfully!');
            } catch (err) {
                console.error("Error importing JagFX:", err);
                alert(`Error importing JagFX: ${err.message}`);
            }
        };
        reader.readAsArrayBuffer(file);
    }

    saveJagFX() {
        this.updateFromUI(); // Ensure model is up-to-date
        const data = this.currentSoundEffect.encode();
        const blob = new Blob([data], { type: 'application/octet-stream' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'sound_effect.jfx';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }

    exportToWav() {
        this.updateFromUI();
        const mixedSamples8Bit = this.currentSoundEffect.mix(); // This returns Int8Array
        if (mixedSamples8Bit.length === 0) {
            alert("No audio data to export.");
            return;
        }

        let audioDataForWav;
        let bitDepthForWav = this.outputBitDepth; // Use selected output bit depth

        if (bitDepthForWav === 8) {
            // Web Audio API expects unsigned 8-bit, so convert signed to unsigned
            audioDataForWav = new Uint8Array(mixedSamples8Bit.length);
            for (let i = 0; i < mixedSamples8Bit.length; i++) {
                audioDataForWav[i] = mixedSamples8Bit[i] + 128;
            }
        } else { // 16-bit
            audioDataForWav = new Int16Array(mixedSamples8Bit.length);
            for (let i = 0; i < mixedSamples8Bit.length; i++) {
                audioDataForWav[i] = mixedSamples8Bit[i] * 256; // Scale 8-bit signed to 16-bit signed
            }
        }

        const wavBuffer = this.encodeWAV(audioDataForWav, 22050, bitDepthForWav);
        const blob = new Blob([wavBuffer], { type: 'audio/wav' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'exported_sound.wav';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }

    // Helper to encode WAV
    encodeWAV(samples, sampleRate, bitDepth) {
        const numChannels = 1;
        const bytesPerSample = bitDepth / 8;
        const blockAlign = numChannels * bytesPerSample;
        const byteRate = sampleRate * blockAlign;
        const dataSize = samples.length * bytesPerSample;

        const buffer = new ArrayBuffer(44 + dataSize);
        const view = new DataView(buffer);

        /* RIFF identifier */
        this._writeString(view, 0, 'RIFF');
        /* RIFF chunk length */
        view.setUint32(4, 36 + dataSize, true);
        /* RIFF type */
        this._writeString(view, 8, 'WAVE');
        /* format chunk identifier */
        this._writeString(view, 12, 'fmt ');
        /* format chunk length */
        view.setUint32(16, 16, true);
        /* sample format (raw) */
        view.setUint16(20, 1, true); // PCM
        /* channel count */
        view.setUint16(22, numChannels, true);
        /* sample rate */
        view.setUint32(24, sampleRate, true);
        /* byte rate (sample rate * block align) */
        view.setUint32(28, byteRate, true);
        /* block align (channel count * bytes per sample) */
        view.setUint16(32, blockAlign, true);
        /* bits per sample */
        view.setUint16(34, bitDepth, true);
        /* data chunk identifier */
        this._writeString(view, 36, 'data');
        /* data chunk length */
        view.setUint32(40, dataSize, true);

        // Write samples
        let offset = 44;
        if (bitDepth === 8) { // samples should be Uint8Array (0-255)
            for (let i = 0; i < samples.length; i++, offset++) {
                view.setUint8(offset, samples[i]);
            }
        } else { // samples should be Int16Array (-32768 to 32767)
            for (let i = 0; i < samples.length; i++, offset += 2) {
                view.setInt16(offset, samples[i], true);
            }
        }
        return buffer;
    }
     _writeString(view, offset, string) {
        for (let i = 0; i < string.length; i++) {
            view.setUint8(offset + i, string.charCodeAt(i));
        }
    }


    playFullEffect() {
        this.updateFromUI();
        const mixedSamples8Bit = this.currentSoundEffect.mix(); // Int8Array
        if (mixedSamples8Bit.length === 0) return;
        this._playAudioData(mixedSamples8Bit, this.outputBitDepth);
    }

    playSelectedTone() {
        this.updateFromUI();
        const tone = this.getSelectedTone();
        if (!tone) return;

        const durationSamples = Math.floor(tone.duration * 22050 / 1000);
        if (durationSamples <= 0) return;

        const synthesized16Bit = tone.synthesize(durationSamples, tone.duration); // Int32Array (effectively 16-bit)

        // Convert to 8-bit signed for consistent playback logic with _playAudioData
        const samples8Bit = new Int8Array(synthesized16Bit.length);
        for(let i=0; i < synthesized16Bit.length; i++) {
            let s = synthesized16Bit[i] >> 8; // Take high byte
            if (s < -128) s = -128; else if (s > 127) s = 127;
            samples8Bit[i] = s;
        }
        this._playAudioData(samples8Bit, this.outputBitDepth);
    }

    _playAudioData(samples8BitSigned, targetPlayBitDepth) {
        this.stopSound();
        if (samples8BitSigned.length === 0) return;

        const audioBuffer = this.audioContext.createBuffer(
            1, // num channels
            samples8BitSigned.length,
            22050 // sample rate
        );
        const channelData = audioBuffer.getChannelData(0); // Float32Array from -1.0 to 1.0

        if (targetPlayBitDepth === 8) {
            for (let i = 0; i < samples8BitSigned.length; i++) {
                channelData[i] = samples8BitSigned[i] / 128.0; // Convert signed 8-bit to float
            }
        } else { // 16-bit
            for (let i = 0; i < samples8BitSigned.length; i++) {
                // Scale the 8-bit sample to represent a 16-bit range then normalize to float
                channelData[i] = (samples8BitSigned[i] * 256) / 32768.0;
            }
        }

        this.currentAudioBufferSource = this.audioContext.createBufferSource();
        this.currentAudioBufferSource.buffer = audioBuffer;
        this.currentAudioBufferSource.connect(this.audioContext.destination);
        this.currentAudioBufferSource.start();
    }


    stopSound() {
        if (this.currentAudioBufferSource) {
            this.currentAudioBufferSource.stop();
            this.currentAudioBufferSource.disconnect();
            this.currentAudioBufferSource = null;
        }
    }

    drawVisualizations() {
        this.drawWaveform();
        const tone = this.getSelectedTone();
        if (tone) {
            this.drawEnvelope(this.pitchEnvCtx, tone.pitch, "Pitch Env");
            this.drawEnvelope(this.volumeEnvCtx, tone.volume, "Volume Env");
            this.drawFilterResponse(this.filterCtx, tone.filter);
        } else {
            this.clearCanvas(this.pitchEnvCtx);
            this.clearCanvas(this.volumeEnvCtx);
            this.clearCanvas(this.filterCtx);
        }
    }

    clearCanvas(ctx) {
        ctx.clearRect(0, 0, ctx.canvas.width, ctx.canvas.height);
        ctx.fillStyle = '#f9f9f9';
        ctx.fillRect(0,0, ctx.canvas.width, ctx.canvas.height);
        ctx.strokeStyle = '#ccc';
        ctx.strokeRect(0,0, ctx.canvas.width, ctx.canvas.height);
    }

    drawWaveform() {
        this.clearCanvas(this.waveformCtx);
        const tone = this.getSelectedTone(); // Or mix full effect for a different view
        if (!tone) return;

        const durationSamples = Math.floor(tone.duration * 22050 / 1000);
        if (durationSamples <= 0) return;

        // Synthesize with current parameters for visualization
        const samples16bit = tone.synthesize(durationSamples, tone.duration);

        const w = this.ui.waveformCanvas.width;
        const h = this.ui.waveformCanvas.height;
        const h2 = h / 2;
        this.waveformCtx.beginPath();
        this.waveformCtx.moveTo(0, h2);
        this.waveformCtx.strokeStyle = '#007bff';

        for (let i = 0; i < durationSamples; i++) {
            const x = (i / durationSamples) * w;
            let yVal;
            if(this.outputBitDepth === 8) {
                yVal = (samples16bit[i] >> 8) / 128.0; // Normalize 8-bit signed
            } else {
                yVal = samples16bit[i] / 32768.0; // Normalize 16-bit signed
            }
            const y = h2 - yVal * h2;
            this.waveformCtx.lineTo(x, y);
        }
        this.waveformCtx.stroke();
    }

    drawEnvelope(ctx, envelope, name) {
        this.clearCanvas(ctx);
        if (!envelope || envelope.form === 0) return;

        const w = ctx.canvas.width;
        const h = ctx.canvas.height;
        ctx.strokeStyle = '#28a745';
        ctx.beginPath();

        const tempEnvelope = new SoundEnvelopeJS(); // Use a temporary copy for simulation
        Object.assign(tempEnvelope, envelope); // Shallow copy properties
        tempEnvelope.durations = [...envelope.durations]; // Deep copy arrays
        tempEnvelope.phases = [...envelope.phases];
        tempEnvelope.reset();

        const graphSteps = 200; // Number of points to plot
        const period = this.getSelectedTone()?.duration || 500; // Use current tone's duration or default

        for (let i = 0; i <= graphSteps; i++) {
            const x = (i / graphSteps) * w;
            const envelopeValue = tempEnvelope.doStep(period); // Simulate step
            // Normalize envelopeValue (0 to 32767 from doStep) to 0-1 range
            const normalizedValue = Math.max(0, Math.min(1, envelopeValue / 32767.0));
            const y = h * (1.0 - normalizedValue); // Invert Y for canvas
            if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
        }
        ctx.stroke();
        ctx.fillStyle = "black";
        ctx.fillText(name, 5, 15);
    }

    drawFilterResponse(ctx, filter) {
        this.clearCanvas(ctx);
         if (!filter || (filter.pairs[0] === 0 && filter.pairs[1] === 0)) {
            ctx.fillStyle = "black";
            ctx.fillText("Filter Off/Flat", 5, 15);
            return;
         }

        const w = ctx.canvas.width;
        const h = ctx.canvas.height;
        ctx.strokeStyle = '#dc3545';
        ctx.beginPath();

        // Simplified: draw a line based on unity values
        // A real filter response graph is much more complex
        const unity0Norm = filter.unity[0] / 65535.0;
        const unity1Norm = filter.unity[1] / 65535.0;

        ctx.moveTo(0, h * (1 - unity0Norm));
        ctx.lineTo(w, h * (1 - unity1Norm));
        ctx.stroke();
        ctx.fillStyle = "black";
        ctx.fillText("Filter (Simplified Unity)", 5, 15);
    }
}

// --- Initialize the Editor ---
document.addEventListener('DOMContentLoaded', () => {
    window.soundEditor = new SoundEffectEditorJS();
});