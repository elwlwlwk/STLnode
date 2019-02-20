class STL{
    constructor(ts, freq, s_window, {s_degree=0, t_window=null,
    t_degree=1, l_window=null, l_degree=null, s_jump=null,
    t_jump=null, l_jump=null, robust=False, inner=null, outer=null}) {

        if(ts[0].length > 2) throw 'The time series must have 1 timestamp and 1 value';
        if(ts.map( datum => datum[1] ).reduce( (a,b) => a+b ) === NaN) throw 'The time series contains NaNs.';

        var n = ts.length;

        if(freq < 2) throw 'The frequency must be greater than 1.';
        if(n <= 2 * freq) throw 'The time series must contain more than 2 full periods of data.';

        if(s_window === 'periodic') s_window = 10 * n + 1;
        if(s_jump === null) s_jump = parseInt(Math.ceil(s_window / 10));

        if(t_window === null) t_window = this.nextodd(parseInt(Math.ceil(1.5 * freq / (1 - 1.5 / s_window))));
        if(t_jump === null) t_jump = parseInt(Math.ceil(t_window / 10));

        if(l_window === null) l_window = this.nextodd(freq);
        if(l_degree === null) l_degree = t_degree;
        if(l_jump === null) l_jump = parseInt(Math.ceil(l_window / 10));

        if(inner === null) robust ? inner = 1 : inner = 2;
        if(outer === null) robust ? outer = 15 : outer = 0;

        var weights = new Array(n).fill(0);
        var seasonal = new Array(n).fill(0);
        var trend = new Array(n).fill(0);
        var work =  new Array(n + 2 * freq).fill(0).map( x => new Array(5).fill(0));

        s_window = Math.max(3, s_window);
        t_window = Math.max(3, t_window);
        l_window = Math.max(3, l_window);
        if(s_window % 2 === 0) s_window += 1;
        if(t_window % 2 === 0) t_window += 1;
        if(l_window % 2 === 0) l_window += 1;

        var userw = false;
        work = this.stlstp(ts, n, freq, s_window, t_window, l_window, s_degree, t_degree, l_degree, s_jump, t_jump, l_jump, inner, userw, weights, seasonal, trend, work);

        userw = true
        for(let _=0; _<outer; _++) {
            var t_ = this.transpose([trend, seasonal]).map( x => x[0]+x[1] );
            for(let i=0; i<n; i++) work[i][0] = t_[i];
            this.stlrwt(ts, n, work.map( x => x[0] ), weights);
            work = this.stlstp(ts, n, freq, s_window, t_window, l_window, s_degree, t_degree, l_degree, s_jump, t_jump, l_jump, inner, userw, weights, seasonal, trend, work);
        }

        if(outer <= 0) weights.fill(1);

        this.seasonal = seasonal
        this.trend = trend
        this.remainder = ts - trend - seasonal
        this.weights = weights

        this.s_window = s_window
        this.t_window = t_window
        this.l_window = l_window
        this.s_degree = s_degree
        this.t_degree = t_degree
        this.l_degree = l_degree
        this.s_jump = s_jump
        this.t_jump = t_jump
        this.l_jump = l_jump
        this.inner = inner
        this.outer = outer
    };

    transpose(m) { return m[0].map((x,i) => m.map(x => x[i])) };

    nextodd(x) {
        x = parseInt(Math.round(x));
        if(x % 2 === 0) x += 1;
        return x;
    };

    stless(y, n, length, ideg, njump, userw, rw, ys, res) {
        if(n < 2) {
            ys[0] = y[0];
            return;
        }

        var newnj = Math.min(njump, n-1);
        if(length >= n) {
            var nleft = 1;
            var nright = n;
            for(let i=0; i<n; i+=newnj) {
                var nys = this.stlest(y, n, length, ideg, i+1, ys[i], nleft, nright, res, userw, rw);
                if(nys !== null) ys[i] = nys;
                else ys[i] = y[i];
            }
        } else {
            if(newnj === 1) {
                var nhs = parseInt((length+1)/2)
                var nleft = 1;
                var nright = length;
                for(let i=0; i<n; i++) {
                    if(i+1 > nhs && nright !== n) {
                        nleft += 1;
                        nright += 1;
                    }
                    let nys = this.stlest(y, n, length, ideg, i+1, ys[i], nleft, nright, res, userw, rw);
                    if(nys !== null) ys[i] = nys;
                    else ys[i] = y[i];
                }
            } else {
                var nhs = parseInt((length+1)/2);
                for(let i=1; i<n+1; i+=newnj) {
                    if(i < nhs) {
                        var nleft = 1;
                        var nright = length;
                    } else if(i >= (n-nhs+1)) {
                        var nleft = n-length+1;
                        var nright = n;
                    } else {
                        var nleft = i-nhs+1;
                        var nright = length+i-nhs;
                    }
                    nys = this.stlest(y, n, length, ideg, i, ys[i-1], nleft, nright, res, userw, rw);
                    if(nys !== null) ys[i-1] = nys;
                    else ys[i-1] = y[i-1];
                }
            }
        }

        if(newnj !== 1) {
            for(let i=0; i<n-newnj; i+=newnj) {
                var delta = (ys[i+newnj] - ys[i]) / newnj;
                for(let j=i+1; j<i+newnj; j++) ys[j] = ys[i] + delta * (j-i);
            }
            var k = parseInt(Math.floor((n-1) / newnj) * newnj + 1);

            if(k !== n) {
                var nys = this.stlest(y, n, length, ideg, n, ys[n-1], nleft, nright, res, userw, rw);
                if(nys !== null) ys[n-1] = nys;
                else ys[n-1] = y[n-1];

                if(k !== n-1) {
                    delta = (ys[n-1] - ys[k-1]) / (n-k);
                    for(let i=k; i<n-1; i++) ys[i] = ys[k-1] + delta * (i - k + 1)
                }
            }
        }
    }

    stlest(y, n, length, ideg, xs, ys, nleft, nright, w, userw, rw) {
        var nleft = parseInt(nleft);
        var nright = parseInt(nright);

        var h = Math.max(xs - nleft, nright - xs);
        if(length > n) h += Math.floor((length - n) / 2);

        var r = [];
        for(let i=nleft-xs; i<nright-xs+1; i++) r.push(Math.abs(i));
        var my_window = [];
        for(let i=nleft-1; i<nright; i++) my_window.push(i);

        var low_mask = r.map( x => x <= 0.001*h );
        var high_mask = r.map( x => x > 0.999*h );
        var mid_mask = this.transpose([low_mask, high_mask]).map( x => !(x[0] || x[1]) );
        var lowmid_mask = high_mask.map( x => !x );

        var low = this.transpose([my_window, low_mask]).filter( x => x[1] ).map( x => x[0] );
        var high = this.transpose([my_window, high_mask]).filter( x => x[1] ).map( x => x[0] );
        var mid = this.transpose([my_window, mid_mask]).filter( x => x[1] ).map( x => x[0] );
        var lowmid = this.transpose([my_window, lowmid_mask]).filter( x => x[1] ).map( x => x[0] );

        var r_mid = this.transpose([r, mid_mask]).filter( x => x[1] ).map( x => x[0] );

        for(let i of low) w[i] = 1;
        for(let i=0; i<mid.length; i++) w[mid[i]] = Math.pow(1 - Math.pow(r_mid[i]/h, 3), 3);
        if(userw) for(let i of lowmid) w[i] *= rw[i];
        var a = lowmid.map( (x) => w[x] ).reduce( (a,b) => a+b );

        for(let i of high) w[i] = 0;

        if(a <= 0) var ret = null;
        else {
            for(let i=nleft-1; i<nright; i++) w[i] /= a;
            if(h > 0 && ideg > 0) {
                a = this.transpose([w.slice(nleft-1, nright), Array.from(Array(nright+1-nleft).keys()).map( (x) => x+nleft )]).map( (x) => x[0] * x[1] ).reduce( (a,b) => a + b );
                var b = xs - a;
                var c = this.transpose([w.slice(nleft-1, nright), Array.from(Array(nright+1-nleft).keys()).map( (x) => Math.pow(x+nleft-a, 2) )]).map( (x) => x[0] * x[1] ).reduce( (a,b) => a + b );
                if(Math.sqrt(c) > 0.001*(n-1)) {
                    b /= c;
                    let t = Array.from(Array(nright+1-nleft).keys()).map( x => (x+nleft-a)*b + 1 );//(b * numpy.arange(nleft-a, nright+1-a) + 1)
                    for(let i=nleft-1; i<nright; i++) w[i] *= t[i-nleft+1];
                }
            }
            ret = this.transpose([w.slice(nleft-1, nright),y.slice(nleft-1, nright)]).map( x => x[0] * x[1] ).reduce( (a,b) => a+b );
        }

        return ret;
    }

    stlfts(x, n, np, trend, work) {
        this.stlma(x, n, np, trend);
        this.stlma(trend, n-np+1, np, work);
        this.stlma(work, n-2*np+2, 3, trend);
    }

    stlma(x, n, length, ave) {
        var v = x.slice(0, length).reduce( (a,b) => a+b );
        ave[0] = v / length;

        var newn = n - length + 1;
        if(newn > 1) {
            var k = length;
            var m = 0;
            for(let j=1; j<newn; j++) {
                k += 1;
                m += 1;
                v = v - x[m-1] + x[k-1];
                ave[j] = v / length;
            }
        }
    }

    stlstp(y, n, np, ns, nt, nl, isdeg, itdeg, ildeg, nsjump, ntjump, nljump, ni, userw, rw, season, trend, work) {
        for(let _=0; _<ni; _++) {
            for(let i=0; i<n; i++) work[i][0] = y[i][1] - trend[i];
            var [work0, work1, work2, work3, work4] = this.transpose(work);
            this.stlss(work0, n, np, ns, isdeg, nsjump, userw, rw, work1, work2, work3, work4, season);
            this.stlfts(work1, n+2*np, np, work2, work0);
            this.stless(work2, n, nl, ildeg, nljump, false, work3, work0, work4);
            for(let i=0; i<n; i++) season[i] = work1[np+i] - work0[i];
            for(let i=0; i<n; i++) work0[i] = y[i][1] - season[i];
            this.stless(work0, n, nt, itdeg, ntjump, userw, rw, trend, work2);
            work = this.transpose([work0, work1, work2, work3, work4]);
        }
        return work;
    }

    stlrwt(y, n, fit, rw) {
        var r = this.transpose([y.map( (x) => x[1]), fit]).map( x => Math.abs(x[0]-x[1]) );
        var sorted_r = r.slice(0).sort( (a,b) => a-b );
        var med = 6 * sorted_r[parseInt(sorted_r.length/2)];
        var low = r.map( x => x <= 0.001*med );
        var high = r.map( x => x > 0.999*med );
        var mid = this.transpose([low, high]).map( x => !(x[0] || x[1]) );

        for(let i in low) low[i] ? rw[i] = 1 : null;
        for(let i in mid) mid[i] ? rw[i] = Math.pow(1 - Math.pow(r[i] / med, 2), 2) : null;
        for(let i in high) high[i] ? rw[i] = 0 : null;
    }

    stlss(y, n, np, ns, isdeg, nsjump, userw, rw, season, work1, work2, work3, work4) {
        for(let j=0; j<np; j++) {
            var k = Math.floor((n-j-1)/np+1);
            for(let i=0; i<k; i++) work1[i] = y[i*np+j];

            if(userw) for(let i=0; i<k; i++) work3[i] = rw[i*np+j];

            var work2_1 = work2.slice(1);
            this.stless(work1, k, ns, isdeg, nsjump, userw, work3, work2_1, work4);
            for(let i=0; i<work2_1.length; i++) work2[i+1] = work2_1[i];
            var nright = Math.min(ns, k);

            var nval = this.stlest(work1, k, ns, isdeg, 0, work2[0], 1, nright, work4, userw, work3);
            if(nval !== null) work2[0] = nval;
            else work2[0] = work2[1];

            var nleft = Math.max(1, k-ns+1);

            nval = this.stlest(work1, k, ns, isdeg, k+1, work2[k+1], nleft, k, work4, userw, work3);
            if(nval !== null) work2[k+1] = nval;
            else work2[k+1] = work2[k];

            for(let m=0; m<k+2; m++) season[m*np+j] = work2[m];
        }
    }
}

module.exports.STL = STL;