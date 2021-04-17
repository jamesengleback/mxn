4: [function(e, t, n) {
    "use strict";
    var r = e("jquery");
    t.exports = function(e, t, n) {
	r.event.addProp("dataTransfer"),
	e.appVersion = "1.13.0",
	e.state = {},
	e.state.input = {},
	e.state.result = {},
	e.state.results1 = [],
	e.state.results2 = [],
	e.state.input.p1 = "",
	e.state.input.p2 = "",
	e.state.input.ct = .25,
	e.state.input.groupKeys = [],
	e.state.input.products = {},
	e.state.input.prodKeys = [],
	e.state.input.group = "",
	e.state.input.product = "",
	e.state.result.tm1 = "---",
	e.state.result.tm2 = "---",
	e.state.result.len1 = "---",
	e.state.result.len2 = "---",
	e.state.result.gc1 = "---",
	e.state.result.gc2 = "---",
	e.state.result.ta = "---",
	e.state.result.itemlist = [],
	e.state.result.critlist = [],
	e.state.p1status = "",
	e.state.p2status = "",
	e.state.results1.showexpanded1 = !1,
	e.state.results2.showexpanded1 = !1;
	var i = window
	  , a = document
	  , o = a.documentElement
	  , s = a.getElementsByTagName("body")[0]
	  , l = i.innerWidth || o.clientWidth || s.clientWidth;
	i.innerHeight || o.clientHeight || s.clientHeight;
	e.initScreenWidth = l,
	e.screenMessage = !1,
	e.openmodal = function(e, n) {
	    var r = {
		templateUrl: "views/modals/" + e + ".html",
		size: n
	    };
	    return t.open(r, {})
	}
	,
	e.$watch(function() {
	    return MathJax.Hub.Queue(["Typeset", MathJax.Hub]),
	    !0
	})
    }
}


    6: [function(e, t, n) {
        "use strict";
        e("angular");
        var r = e("../vendorlibs/FileSaver");
        t.exports = function(e, t, n, i, a, o, s, l, u) {
            e.input = {},
            e.result = {},
            e.output = [],
            e.input.p1 = "",
            e.input.p2 = "",
            e.input.id1 = "",
            e.input.id2 = "",
            e.input.ct = .25,
            e.result.tm1 = "",
            e.result.tm2 = "",
            e.result.ta = "",
            e.result.itemlist = [],
            e.result.critlist = [],
            e.input.groupKeys = [],
            e.input.products = {},
            e.input.prodKeys = [],
            e.input.group = "",
            e.input.product = "",
            e.input.batch = "",
            e.input.interleaved = !1,
            e.input.filename = "",
            e.result.batch = "",
            e.result.batch2 = "",
            e.p1status = "",
            e.p2status = "",
            console.log(n),
            i.setData(n),
            console.log(i.tmcdata),
            i.restoreUserPrefs() ? (e.lastp1 = i.restoreUserPrefs().p1,
            e.lastp2 = i.restoreUserPrefs().p2) : (e.lastp1 = "",
            e.lastp2 = "");
            var c = s.p1seq || e.lastp1 || ""
              , p = s.p2seq || e.lastp2 || "";
            e.prefill = function() {
                e.input.interleaved ? e.input.batch = "P1fwd; AGCGGATAACAATTTCACACAGGA\nP1rev; GTAAAACGACGGCCAGT\nP3fwd; AGCGGATAAGGGCAATTTCAC\nP3rev; GTAAAACGACGGCCA\n" : e.input.batch = "P1fwd; AGCGGATAACAATTTCACACAGGA; P1rev; GTAAAACGACGGCCAGT\npBeta; AGCGGATAACAATTTCAC\nP3fwd; AGCGGATAAGGGCAATTTCAC; P3rev; GTAAAACGACGGCCA\n",
                e.output.showresultstable = !1,
                e.data_toggle_label = e.output.showresultstable ? "Hide" : "Show",
                e.runCalc3(e.input.interleaved)
            }
            ,
            e.clearCalc = function() {
                e.input.batch = "",
                e.input.filename = "",
                e.output = [],
                e.result.itemlist = [],
                e.result.critlist = [],
                e.result.batch = "",
                e.result.batch2 = "",
                e.runmsg = "",
                document.getElementById("fileinput").value = "",
                e.output.showresultstable = !1,
                e.data_toggle_label = e.output.showresultstable ? "Hide" : "Show"
            }
            ,
            e.setGroups = function() {
                e.input.groupKeys = i.getGroupKeys(),
                e.input.group = e.input.groupKeys[0]
            }
            ,
            e.setProducts = function() {
                e.input.group;
                e.input.productKeys = i.getProductKeysForGroup(e.input.group),
                e.input.products = i.getProductsForGroup(e.input.group),
                e.input.product = e.input.products[0].id;
                e.input.product;
                e.setCt()
            }
            ,
            e.setCt = function() {
                var t = e.input.group
                  , n = e.input.product;
                0 === t.indexOf("Phusion") || 0 === t.indexOf("Q5") ? e.input.ct = 500 : 0 === t.indexOf("LongAmp") ? e.input.ct = 400 : 0 === t.indexOf("Master") ? 0 === n.indexOf("phusion") || 0 === n.indexOf("q5") ? e.input.ct = 500 : 0 === n.indexOf("lataq") || 0 === n.indexOf("lahstaq") ? e.input.ct = 400 : e.input.ct = 200 : e.input.ct = 200,
                e.runCalc3(e.input.interleaved)
            }
            ,
            e.runCalc3 = function(t) {
                t = !!t,
                console.log("Interleaved", t);
                for (var n, r, s, l, u, c, p, d, f, h, g, m, v, $, b, y, w, x, C, k = "", S = "", T = "", A = "", D = "", M = e.input.ct / 1e3, E = e.input.product, O = e.input.group, N = [], P = [], q = 0, I = null, R = 1, L = e.input.batch, V = L.split(/\n\r?/), j = 0; j < V.length; ++j)
                    V[j].match(/^\s*$/) ? V.splice(j, 1) : V[j] = V[j].replace(/[;,\t\s]$/, "");
                if (V && 0 !== V.length) {
                    console.log("after cleanup");
                    for (var _ = V[0].slice(), F = _.split(""), H = 0; H < F.length; ++H)
                        console.log(H, F[H], _.charCodeAt(H));
                    var U = {};
                    if (console.log("testing separators"),
                    U.semicolon = _.split(";").length,
                    console.log("fields when split with semicolon: ", U.semicolon),
                    U.comma = _.split(",").length,
                    console.log("fields when split with comma: ", U.comma),
                    U.tab = _.split("\t").length,
                    console.log("fields when split with tab: ", U.tab),
                    U.ws = _.split(/[\s\W]+/).length,
                    console.log("fields when split with ws: ", U.ws),
                    t ? (I = 2 === U.semicolon ? ";" : 2 === U.comma ? "," : 2 === U.tab ? "\t" : 2 === U.ws ? /[\s\W]+/ : null,
                    R = 2) : (I = 2 === U.semicolon || 3 === U.semicolon || 4 === U.semicolon ? ";" : 2 === U.comma || 3 === U.comma || 4 === U.comma ? "," : 2 === U.tab || 3 === U.tab || 4 === U.tab ? "\t" : 2 === U.ws || 3 === U.ws || 4 === U.ws ? /[\s\W]+/ : null,
                    R = 1),
                    null == I)
                        return P = ["Unable to read primers and/or invalid data format. "],
                        e.result.critlist = P,
                        e.result.itemlist = [],
                        L = "",
                        e.input.filename = "",
                        document.getElementById("fileinput").value = "",
                        e.runmsg = "",
                        void (e.output = []);
                    if (e.result = {},
                    e.output = [],
                    e.ctstatus = "",
                    e.runmsg = "",
                    e.output.showresultstable = !1,
                    e.data_toggle_label = e.output.showresultstable ? "Hide" : "Show",
                    M <= 0 || isNaN(M / 1))
                        return !1,
                        e.ctstatus = "invalidct",
                        P.push("Invalid primer concentration. "),
                        e.result.itemlist = N,
                        e.result.critlist = P,
                        void (e.p1status = "");
                    switch ("onetaq_gc" === (s = i.getBufferIdForProduct(E)) && (q = 5),
                    l = i.getBufferSaltForProduct(e.input.product),
                    O) {
                    case "Phusion":
                    case "Phusion Hot Start Flex":
                        u = 5,
                        M /= 4;
                        break;
                    case "Master Mix":
                        0 === E.indexOf("phusion") ? (u = 5,
                        M /= 4) : u = 4;
                        break;
                    default:
                        u = 4
                    }
                    ($ = i.validateBuffer(E, O, s)).hasCritWarnings && Array.prototype.push.apply(P, $.critwarnings),
                    $.hasWarnings && Array.prototype.push.apply(N, $.warnings),
                    e.result.itemlist = N,
                    e.result.critlist = P,
                    a.setCt(M).setMonosalt(l).setMethod(u).setDMSO(q),
                    p = 0,
                    d = 0,
                    g = 0;
                    for (var B = 0; B <= V.length - R; B += R) {
                        if (D = "OK",
                        f = !1,
                        h = !1,
                        t) {
                            if (n = V[B].split(I),
                            r = V[B + 1].split(I),
                            n.length < 2 && r.length < 2) {
                                ++g;
                                continue
                            }
                            2 === n.length ? (k = n[1],
                            T = n[0]) : (k = "",
                            T = ""),
                            2 === r.length ? (S = r[1],
                            A = r[0]) : (S = "",
                            A = "")
                        } else if (4 === (n = V[B].split(I)).length)
                            k = n[1].replace(/[\/W\/s]/g, ""),
                            S = n[3].replace(/[\/W\/s]/g, ""),
                            T = n[0],
                            A = n[2];
                        else if (3 === n.length)
                            k = n[1].replace(/[\/W\/s]/g, ""),
                            S = n[2].replace(/[\/W\/s]/g, ""),
                            T = n[0],
                            A = n[0];
                        else {
                            if (2 !== n.length) {
                                ++g;
                                continue
                            }
                            k = n[1].replace(/[\/W\/s]/g, ""),
                            S = "",
                            T = n[0],
                            A = ""
                        }
                        if (e.input.p1 = k.toUpperCase(),
                        e.input.p2 = S.toUpperCase(),
                        e.input.id1 = T,
                        e.input.id2 = A,
                        (m = i.validateInput(e.input.p1.replace(/\s/g, ""), e.input.p2.replace(/\s/g, ""), M, E, O)).p1isValid,
                        m.p2isValid,
                        c = m.hasCritWarnings,
                        e.p1status = m.p1status,
                        e.p2status = m.p2status,
                        c && (D = m.critwarnings.join("-- "),
                        f = !0),
                        "invalidseq" !== e.p1status) {
                            var G, W, z = o.expandDNASeq(k.replace(/\s/g, ""));
                            if (1 === z.length)
                                y = (W = a.setCt(M).setMonosalt(l).setMethod(u).setDMSO(q).calcTm_nomm(k.replace(/\s/g, ""))).tm,
                                y = Math.round(Math.round(10 * y) / 10),
                                e.result.tm1 = y,
                                e.result.len1 = W.len,
                                e.result.gc1 = 100 * W.fgc,
                                e.result.gc1 = Math.round(10 * e.result.gc1 / 10);
                            else {
                                e.results1 = [],
                                a.setCt(M).setMonosalt(l).setMethod(u).setDMSO(q);
                                for (var Y = 0; Y < z.length; ++Y)
                                    G = z[Y],
                                    W = a.calcTm_nomm(G),
                                    y = Math.round(Math.round(10 * W.tm) / 10),
                                    e.results1.push({
                                        seq: G,
                                        tm: y,
                                        len: G.length,
                                        gc: Math.round(100 * W.fgc * 10 / 10)
                                    });
                                y = e.results1.reduce(function(e, t) {
                                    return t.tm < e && (e = t.tm),
                                    e
                                }, e.results1[0].tm),
                                x = e.results1.reduce(function(e, t) {
                                    return t.tm > e && (e = t.tm),
                                    e
                                }, e.results1[0].tm),
                                e.result.tm1 = y !== x ? y + " - " + x : y,
                                e.result.len1 = k.replace(/\s/g, "").length,
                                e.result.gc1 = "---",
                                e.result.tm1m = x
                            }
                        }
                        if ("invalidseq" !== e.p2status) {
                            var K, Q, X = o.expandDNASeq(S.replace(/\s/g, ""));
                            if (1 === X.length)
                                w = (Q = a.setCt(M).setMonosalt(l).setMethod(u).setDMSO(q).calcTm_nomm(S.replace(/\s/g, ""))).tm,
                                w = Math.round(Math.round(10 * w) / 10),
                                e.result.tm2 = w,
                                e.result.len2 = Q.len,
                                e.result.gc2 = 100 * Q.fgc,
                                e.result.gc2 = Math.round(10 * e.result.gc2 / 10);
                            else {
                                e.results2 = [],
                                a.setCt(M).setMonosalt(l).setMethod(u).setDMSO(q);
                                for (var J = 0; J < z.length; ++J)
                                    K = X[J],
                                    Q = a.calcTm_nomm(K),
                                    w = Math.round(Math.round(10 * Q.tm) / 10),
                                    e.results2.push({
                                        seq: K,
                                        tm: w,
                                        len: K.length,
                                        gc: Math.round(100 * Q.fgc * 10 / 10)
                                    });
                                w = e.results2.reduce(function(e, t) {
                                    return t.tm < e && (e = t.tm),
                                    e
                                }, e.results2[0].tm),
                                C = e.results2.reduce(function(e, t) {
                                    return t.tm > e && (e = t.tm),
                                    e
                                }, e.results2[0].tm),
                                e.result.tm2 = w !== C ? w + " - " + C : w,
                                e.result.len2 = S.replace(/\s/g, "").length,
                                e.result.gc2 = "---",
                                e.result.tm2m = C
                            }
                        }
                        "invalidseq" !== e.p1status && "invalidseq" !== e.p2status ? (b = i.getAnnealTemp(k, y, S, w, O, E),
                        e.result.ta = Math.round(b),
                        (v = i.validateTm(y, w, b, E, O, s)).hasCritWarnings && (h = !0,
                        D += "--" + v.critwarnings.join("--")),
                        v.hasWarnings && (h = !0,
                        D += "-- " + v.warnings.join("-- "))) : e.result.ta = "---",
                        f && (p += 1),
                        h && (d += 1),
                        "invalidseq" !== e.p1status && e.output.push([e.input.id1, e.input.p1, e.result.tm1, e.result.ta, D].join("\t")),
                        "invalidseq" !== e.p2status && e.output.push([e.input.id2, e.input.p2, e.result.tm2, e.result.ta, D].join("\t")),
                        e.p1status = ""
                    }
                    e.output.length,
                    e.runmsg = g > 0 ? g + " Invalid line(s) " : "",
                    e.runmsg += e.output.length + " primers processed. Errors: " + p + " Warnings: " + d,
                    e.result.batch = e.output.join("\n"),
                    e.result.batch2 = [];
                    for (var Z = 0; Z < e.output.length; ++Z)
                        e.result.batch2.push(e.output[Z].split("\t"));
                    0 == e.output.length && V.length > 0 ? e.novaliddatamsg = "Unable to detect format or no valid data entered -- check format of 1st line." : e.novaliddatamsg = ""
                }
            }
            ,
            e.getBatchResults = function() {
                return e.result.batch
            }
            ,
            e.downloadData = function() {
                var t = e.getBatchResults()
                  , n = new Blob([t],{
                    type: "text/plain;charset=utf-8"
                });
                r.saveAs(n, "tmcalc_batch.txt")
            }
            ;
            e.$on("$routeChangeStart", function() {}),
            e.about = function() {
                l.path("#/about")
            }
            ,
            e.safeApply = function(e) {
                var t = this.$root.$$phase;
                "$apply" == t || "$digest" == t ? e && "function" == typeof e && e() : this.$apply(e)
            }
            ,
            e.processFileData = function(t) {
                e.safeApply(function() {
                    e.input.batch = t,
                    console.log("filedata:" + t),
                    e.output.showresultstable = !1,
                    e.data_toggle_label = e.output.showresultstable ? "Hide" : "Show",
                    e.runCalc3(e.input.interleaved)
                })
            }
            ,
            e.switch2single = function() {
                l.path("/")
            }
            ,
            e.toggle_data_display = function() {
                e.output.showresultstable = !e.output.showresultstable,
                e.data_toggle_label = e.output.showresultstable ? "Hide" : "Show"
            }
            ,
            e.output.showresultstable = !1,
            e.data_toggle_label = e.output.showresultstable ? "Hide" : "Show",
            e.input.p1 = c,
            e.input.p2 = p,
            e.setGroups(),
            e.setProducts()
        }
    }
    , {
        "../vendorlibs/FileSaver": 22,
        angular: 33
    }],


    7: [function(e, t, n) {
        var r = e("angular").module("tmcApp");
        r.controller("MainCtrl", ["$scope", "$rootScope", "nebutil", "prodData", "tmcalculatorData", "tmcalc", "neb_bioseq", "$routeParams", "$location", "$window", e("./main")]),
        r.controller("BatchCtrl", ["$scope", "nebutil", "prodData", "tmcalculatorData", "tmcalc", "neb_bioseq", "$routeParams", "$location", "$window", e("./batch")]),
        r.controller("AboutCtrl", ["$scope", e("./about")])
    }
    , {
        "./about": 5,
        "./batch": 6,
        "./main": 8,
        angular: 33
    }],
