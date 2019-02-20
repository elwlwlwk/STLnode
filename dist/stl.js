'use strict';

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

var STL = function () {
    function STL(ts, freq, s_window, _ref) {
        var _ref$s_degree = _ref.s_degree,
            s_degree = _ref$s_degree === undefined ? 0 : _ref$s_degree,
            _ref$t_window = _ref.t_window,
            t_window = _ref$t_window === undefined ? null : _ref$t_window,
            _ref$t_degree = _ref.t_degree,
            t_degree = _ref$t_degree === undefined ? 1 : _ref$t_degree,
            _ref$l_window = _ref.l_window,
            l_window = _ref$l_window === undefined ? null : _ref$l_window,
            _ref$l_degree = _ref.l_degree,
            l_degree = _ref$l_degree === undefined ? null : _ref$l_degree,
            _ref$s_jump = _ref.s_jump,
            s_jump = _ref$s_jump === undefined ? null : _ref$s_jump,
            _ref$t_jump = _ref.t_jump,
            t_jump = _ref$t_jump === undefined ? null : _ref$t_jump,
            _ref$l_jump = _ref.l_jump,
            l_jump = _ref$l_jump === undefined ? null : _ref$l_jump,
            _ref$robust = _ref.robust,
            robust = _ref$robust === undefined ? False : _ref$robust,
            _ref$inner = _ref.inner,
            inner = _ref$inner === undefined ? null : _ref$inner,
            _ref$outer = _ref.outer,
            outer = _ref$outer === undefined ? null : _ref$outer;

        _classCallCheck(this, STL);

        if (ts[0].length > 2) throw 'The time series must have 1 timestamp and 1 value';
        if (ts.map(function (datum) {
            return datum[1];
        }).reduce(function (a, b) {
            return a + b;
        }) === NaN) throw 'The time series contains NaNs.';

        var n = ts.length;

        if (freq < 2) throw 'The frequency must be greater than 1.';
        if (n <= 2 * freq) throw 'The time series must contain more than 2 full periods of data.';

        if (s_window === 'periodic') s_window = 10 * n + 1;
        if (s_jump === null) s_jump = parseInt(Math.ceil(s_window / 10));

        if (t_window === null) t_window = this.nextodd(parseInt(Math.ceil(1.5 * freq / (1 - 1.5 / s_window))));
        if (t_jump === null) t_jump = parseInt(Math.ceil(t_window / 10));

        if (l_window === null) l_window = this.nextodd(freq);
        if (l_degree === null) l_degree = t_degree;
        if (l_jump === null) l_jump = parseInt(Math.ceil(l_window / 10));

        if (inner === null) robust ? inner = 1 : inner = 2;
        if (outer === null) robust ? outer = 15 : outer = 0;

        var weights = new Array(n).fill(0);
        var seasonal = new Array(n).fill(0);
        var trend = new Array(n).fill(0);
        var work = new Array(n + 2 * freq).fill(0).map(function (x) {
            return new Array(5).fill(0);
        });

        s_window = Math.max(3, s_window);
        t_window = Math.max(3, t_window);
        l_window = Math.max(3, l_window);
        if (s_window % 2 === 0) s_window += 1;
        if (t_window % 2 === 0) t_window += 1;
        if (l_window % 2 === 0) l_window += 1;

        var userw = false;
        work = this.stlstp(ts, n, freq, s_window, t_window, l_window, s_degree, t_degree, l_degree, s_jump, t_jump, l_jump, inner, userw, weights, seasonal, trend, work);

        userw = true;
        for (var _ = 0; _ < outer; _++) {
            var t_ = this.transpose([trend, seasonal]).map(function (x) {
                return x[0] + x[1];
            });
            for (var i = 0; i < n; i++) {
                work[i][0] = t_[i];
            }this.stlrwt(ts, n, work.map(function (x) {
                return x[0];
            }), weights);
            work = this.stlstp(ts, n, freq, s_window, t_window, l_window, s_degree, t_degree, l_degree, s_jump, t_jump, l_jump, inner, userw, weights, seasonal, trend, work);
        }

        if (outer <= 0) weights.fill(1);

        this.seasonal = seasonal;
        this.trend = trend;
        this.remainder = ts - trend - seasonal;
        this.weights = weights;

        this.s_window = s_window;
        this.t_window = t_window;
        this.l_window = l_window;
        this.s_degree = s_degree;
        this.t_degree = t_degree;
        this.l_degree = l_degree;
        this.s_jump = s_jump;
        this.t_jump = t_jump;
        this.l_jump = l_jump;
        this.inner = inner;
        this.outer = outer;
    }

    _createClass(STL, [{
        key: 'transpose',
        value: function transpose(m) {
            return m[0].map(function (x, i) {
                return m.map(function (x) {
                    return x[i];
                });
            });
        }
    }, {
        key: 'nextodd',
        value: function nextodd(x) {
            x = parseInt(Math.round(x));
            if (x % 2 === 0) x += 1;
            return x;
        }
    }, {
        key: 'stless',
        value: function stless(y, n, length, ideg, njump, userw, rw, ys, res) {
            if (n < 2) {
                ys[0] = y[0];
                return;
            }

            var newnj = Math.min(njump, n - 1);
            if (length >= n) {
                var nleft = 1;
                var nright = n;
                for (var i = 0; i < n; i += newnj) {
                    var nys = this.stlest(y, n, length, ideg, i + 1, ys[i], nleft, nright, res, userw, rw);
                    if (nys !== null) ys[i] = nys;else ys[i] = y[i];
                }
            } else {
                if (newnj === 1) {
                    var nhs = parseInt((length + 1) / 2);
                    var nleft = 1;
                    var nright = length;
                    for (var _i = 0; _i < n; _i++) {
                        if (_i + 1 > nhs && nright !== n) {
                            nleft += 1;
                            nright += 1;
                        }
                        var _nys = this.stlest(y, n, length, ideg, _i + 1, ys[_i], nleft, nright, res, userw, rw);
                        if (_nys !== null) ys[_i] = _nys;else ys[_i] = y[_i];
                    }
                } else {
                    var nhs = parseInt((length + 1) / 2);
                    for (var _i2 = 1; _i2 < n + 1; _i2 += newnj) {
                        if (_i2 < nhs) {
                            var nleft = 1;
                            var nright = length;
                        } else if (_i2 >= n - nhs + 1) {
                            var nleft = n - length + 1;
                            var nright = n;
                        } else {
                            var nleft = _i2 - nhs + 1;
                            var nright = length + _i2 - nhs;
                        }
                        nys = this.stlest(y, n, length, ideg, _i2, ys[_i2 - 1], nleft, nright, res, userw, rw);
                        if (nys !== null) ys[_i2 - 1] = nys;else ys[_i2 - 1] = y[_i2 - 1];
                    }
                }
            }

            if (newnj !== 1) {
                for (var _i3 = 0; _i3 < n - newnj; _i3 += newnj) {
                    var delta = (ys[_i3 + newnj] - ys[_i3]) / newnj;
                    for (var j = _i3 + 1; j < _i3 + newnj; j++) {
                        ys[j] = ys[_i3] + delta * (j - _i3);
                    }
                }
                var k = parseInt(Math.floor((n - 1) / newnj) * newnj + 1);

                if (k !== n) {
                    var nys = this.stlest(y, n, length, ideg, n, ys[n - 1], nleft, nright, res, userw, rw);
                    if (nys !== null) ys[n - 1] = nys;else ys[n - 1] = y[n - 1];

                    if (k !== n - 1) {
                        delta = (ys[n - 1] - ys[k - 1]) / (n - k);
                        for (var _i4 = k; _i4 < n - 1; _i4++) {
                            ys[_i4] = ys[k - 1] + delta * (_i4 - k + 1);
                        }
                    }
                }
            }
        }
    }, {
        key: 'stlest',
        value: function stlest(y, n, length, ideg, xs, ys, nleft, nright, w, userw, rw) {
            var nleft = parseInt(nleft);
            var nright = parseInt(nright);

            var h = Math.max(xs - nleft, nright - xs);
            if (length > n) h += Math.floor((length - n) / 2);

            var r = [];
            for (var i = nleft - xs; i < nright - xs + 1; i++) {
                r.push(Math.abs(i));
            }var my_window = [];
            for (var _i5 = nleft - 1; _i5 < nright; _i5++) {
                my_window.push(_i5);
            }var low_mask = r.map(function (x) {
                return x <= 0.001 * h;
            });
            var high_mask = r.map(function (x) {
                return x > 0.999 * h;
            });
            var mid_mask = this.transpose([low_mask, high_mask]).map(function (x) {
                return !(x[0] || x[1]);
            });
            var lowmid_mask = high_mask.map(function (x) {
                return !x;
            });

            var low = this.transpose([my_window, low_mask]).filter(function (x) {
                return x[1];
            }).map(function (x) {
                return x[0];
            });
            var high = this.transpose([my_window, high_mask]).filter(function (x) {
                return x[1];
            }).map(function (x) {
                return x[0];
            });
            var mid = this.transpose([my_window, mid_mask]).filter(function (x) {
                return x[1];
            }).map(function (x) {
                return x[0];
            });
            var lowmid = this.transpose([my_window, lowmid_mask]).filter(function (x) {
                return x[1];
            }).map(function (x) {
                return x[0];
            });

            var r_mid = this.transpose([r, mid_mask]).filter(function (x) {
                return x[1];
            }).map(function (x) {
                return x[0];
            });

            var _iteratorNormalCompletion = true;
            var _didIteratorError = false;
            var _iteratorError = undefined;

            try {
                for (var _iterator = low[Symbol.iterator](), _step; !(_iteratorNormalCompletion = (_step = _iterator.next()).done); _iteratorNormalCompletion = true) {
                    var _i10 = _step.value;
                    w[_i10] = 1;
                }
            } catch (err) {
                _didIteratorError = true;
                _iteratorError = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion && _iterator.return) {
                        _iterator.return();
                    }
                } finally {
                    if (_didIteratorError) {
                        throw _iteratorError;
                    }
                }
            }

            for (var _i6 = 0; _i6 < mid.length; _i6++) {
                w[mid[_i6]] = Math.pow(1 - Math.pow(r_mid[_i6] / h, 3), 3);
            }if (userw) {
                var _iteratorNormalCompletion2 = true;
                var _didIteratorError2 = false;
                var _iteratorError2 = undefined;

                try {
                    for (var _iterator2 = lowmid[Symbol.iterator](), _step2; !(_iteratorNormalCompletion2 = (_step2 = _iterator2.next()).done); _iteratorNormalCompletion2 = true) {
                        var _i7 = _step2.value;
                        w[_i7] *= rw[_i7];
                    }
                } catch (err) {
                    _didIteratorError2 = true;
                    _iteratorError2 = err;
                } finally {
                    try {
                        if (!_iteratorNormalCompletion2 && _iterator2.return) {
                            _iterator2.return();
                        }
                    } finally {
                        if (_didIteratorError2) {
                            throw _iteratorError2;
                        }
                    }
                }
            }var a = lowmid.map(function (x) {
                return w[x];
            }).reduce(function (a, b) {
                return a + b;
            });

            var _iteratorNormalCompletion3 = true;
            var _didIteratorError3 = false;
            var _iteratorError3 = undefined;

            try {
                for (var _iterator3 = high[Symbol.iterator](), _step3; !(_iteratorNormalCompletion3 = (_step3 = _iterator3.next()).done); _iteratorNormalCompletion3 = true) {
                    var _i11 = _step3.value;
                    w[_i11] = 0;
                }
            } catch (err) {
                _didIteratorError3 = true;
                _iteratorError3 = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion3 && _iterator3.return) {
                        _iterator3.return();
                    }
                } finally {
                    if (_didIteratorError3) {
                        throw _iteratorError3;
                    }
                }
            }

            if (a <= 0) var ret = null;else {
                for (var _i8 = nleft - 1; _i8 < nright; _i8++) {
                    w[_i8] /= a;
                }if (h > 0 && ideg > 0) {
                    a = this.transpose([w.slice(nleft - 1, nright), Array.from(Array(nright + 1 - nleft).keys()).map(function (x) {
                        return x + nleft;
                    })]).map(function (x) {
                        return x[0] * x[1];
                    }).reduce(function (a, b) {
                        return a + b;
                    });
                    var b = xs - a;
                    var c = this.transpose([w.slice(nleft - 1, nright), Array.from(Array(nright + 1 - nleft).keys()).map(function (x) {
                        return Math.pow(x + nleft - a, 2);
                    })]).map(function (x) {
                        return x[0] * x[1];
                    }).reduce(function (a, b) {
                        return a + b;
                    });
                    if (Math.sqrt(c) > 0.001 * (n - 1)) {
                        b /= c;
                        var t = Array.from(Array(nright + 1 - nleft).keys()).map(function (x) {
                            return (x + nleft - a) * b + 1;
                        }); //(b * numpy.arange(nleft-a, nright+1-a) + 1)
                        for (var _i9 = nleft - 1; _i9 < nright; _i9++) {
                            w[_i9] *= t[_i9 - nleft + 1];
                        }
                    }
                }
                ret = this.transpose([w.slice(nleft - 1, nright), y.slice(nleft - 1, nright)]).map(function (x) {
                    return x[0] * x[1];
                }).reduce(function (a, b) {
                    return a + b;
                });
            }

            return ret;
        }
    }, {
        key: 'stlfts',
        value: function stlfts(x, n, np, trend, work) {
            this.stlma(x, n, np, trend);
            this.stlma(trend, n - np + 1, np, work);
            this.stlma(work, n - 2 * np + 2, 3, trend);
        }
    }, {
        key: 'stlma',
        value: function stlma(x, n, length, ave) {
            var v = x.slice(0, length).reduce(function (a, b) {
                return a + b;
            });
            ave[0] = v / length;

            var newn = n - length + 1;
            if (newn > 1) {
                var k = length;
                var m = 0;
                for (var j = 1; j < newn; j++) {
                    k += 1;
                    m += 1;
                    v = v - x[m - 1] + x[k - 1];
                    ave[j] = v / length;
                }
            }
        }
    }, {
        key: 'stlstp',
        value: function stlstp(y, n, np, ns, nt, nl, isdeg, itdeg, ildeg, nsjump, ntjump, nljump, ni, userw, rw, season, trend, work) {
            for (var _ = 0; _ < ni; _++) {
                for (var i = 0; i < n; i++) {
                    work[i][0] = y[i][1] - trend[i];
                }
                var _transpose = this.transpose(work),
                    _transpose2 = _slicedToArray(_transpose, 5),
                    work0 = _transpose2[0],
                    work1 = _transpose2[1],
                    work2 = _transpose2[2],
                    work3 = _transpose2[3],
                    work4 = _transpose2[4];

                this.stlss(work0, n, np, ns, isdeg, nsjump, userw, rw, work1, work2, work3, work4, season);
                this.stlfts(work1, n + 2 * np, np, work2, work0);
                this.stless(work2, n, nl, ildeg, nljump, false, work3, work0, work4);
                for (var _i12 = 0; _i12 < n; _i12++) {
                    season[_i12] = work1[np + _i12] - work0[_i12];
                }for (var _i13 = 0; _i13 < n; _i13++) {
                    work0[_i13] = y[_i13][1] - season[_i13];
                }this.stless(work0, n, nt, itdeg, ntjump, userw, rw, trend, work2);
                work = this.transpose([work0, work1, work2, work3, work4]);
            }
            return work;
        }
    }, {
        key: 'stlrwt',
        value: function stlrwt(y, n, fit, rw) {
            var r = this.transpose([y.map(function (x) {
                return x[1];
            }), fit]).map(function (x) {
                return Math.abs(x[0] - x[1]);
            });
            var sorted_r = r.slice(0).sort(function (a, b) {
                return a - b;
            });
            var med = 6 * sorted_r[parseInt(sorted_r.length / 2)];
            var low = r.map(function (x) {
                return x <= 0.001 * med;
            });
            var high = r.map(function (x) {
                return x > 0.999 * med;
            });
            var mid = this.transpose([low, high]).map(function (x) {
                return !(x[0] || x[1]);
            });

            for (var i in low) {
                low[i] ? rw[i] = 1 : null;
            }for (var _i14 in mid) {
                mid[_i14] ? rw[_i14] = Math.pow(1 - Math.pow(r[_i14] / med, 2), 2) : null;
            }for (var _i15 in high) {
                high[_i15] ? rw[_i15] = 0 : null;
            }
        }
    }, {
        key: 'stlss',
        value: function stlss(y, n, np, ns, isdeg, nsjump, userw, rw, season, work1, work2, work3, work4) {
            for (var j = 0; j < np; j++) {
                var k = Math.floor((n - j - 1) / np + 1);
                for (var i = 0; i < k; i++) {
                    work1[i] = y[i * np + j];
                }if (userw) for (var _i16 = 0; _i16 < k; _i16++) {
                    work3[_i16] = rw[_i16 * np + j];
                }var work2_1 = work2.slice(1);
                this.stless(work1, k, ns, isdeg, nsjump, userw, work3, work2_1, work4);
                for (var _i17 = 0; _i17 < work2_1.length; _i17++) {
                    work2[_i17 + 1] = work2_1[_i17];
                }var nright = Math.min(ns, k);

                var nval = this.stlest(work1, k, ns, isdeg, 0, work2[0], 1, nright, work4, userw, work3);
                if (nval !== null) work2[0] = nval;else work2[0] = work2[1];

                var nleft = Math.max(1, k - ns + 1);

                nval = this.stlest(work1, k, ns, isdeg, k + 1, work2[k + 1], nleft, k, work4, userw, work3);
                if (nval !== null) work2[k + 1] = nval;else work2[k + 1] = work2[k];

                for (var m = 0; m < k + 2; m++) {
                    season[m * np + j] = work2[m];
                }
            }
        }
    }]);

    return STL;
}();

module.exports.STL = STL;
//# sourceMappingURL=stl.js.map