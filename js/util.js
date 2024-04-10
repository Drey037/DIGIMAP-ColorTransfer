// namespace slib
var slib = slib || {};


// namespace slib.internal
slib.internal = slib.internal || {};


// types
slib.AssertionError = function(msg)
{
    this.message = msg;
    this.name    = 'AssertionError';
}


slib.ValueError = function(msg)
{
    this.message = msg;
    this.name    = 'ValueError';
}


slib.IndexError = function(msg)
{
    this.message = msg;
    this.name    = 'IndexError';
}

// functions
// inheritance
Function.prototype.inheritFrom = function(SuperClass)
{
    // from: http://phrogz.net/JS/classes/OOPinJS2.html
    if (SuperClass.constructor === Function) {
        this.prototype = new SuperClass;
        this.prototype.constructor = this;
        this.prototype.SuperClass = SuperClass.prototype;
    } else {
        this.prototype = SuperClass;
        this.prototype.constructor = this;
        this.prototype.SuperClass = SuperClass;
    }
}

// assertions
slib.assert = function(expr, errMsg)
{
    if (!expr) {
        errMsg = errMsg || 'Assertion Failed.';
        throw new slib.AssertionError(errMsg);
    }
}

// string formatting
String.prototype.format = function(/* ... */)
{
    var args = arguments;
    return this.replace(/{(\d+)}/g, function(s, i) {
        if (i < args.length) {
            return args[i];
        } else {
            throw slib.IndexError('Index out of range.');
        }
    });
}

// for-loop
slib.foreach = function(that, arr, func)
{
    for (var i = 0; i < arr.length; i++) {
        func.call(that, arr[i]);
    }
}

slib.internal.appendRange = function(start, stop, step)
{
    if (arguments.length == 1) {
        stop  = arguments[0];
        start = 0;
    }
    if (arguments.length < 3) {
        step  = 1;
    }
    
    if (step > 0) {
        for (var i = start; i < stop; i += step) {
            this.push(i);
        }
    } else {
        for (var i = start; i > stop; i += step) {
            this.push(i);
        }
    }
    return this;
}

slib.range = function(start, stop, step)
{
    return slib.internal.appendRange.apply([], arguments);
}

slib.enumerate = function(arr)
{
    var result = [];
    for (var i = 0; i < arr.length; i++) {
        result.push([i, arr[i]]);
    }
    return result;
}

slib.zip = function(/* ... */)
{
    var result = [];
    if (arguments.length > 0) {
        var n = arguments[0].length;
        for (var i = 1; i < arguments.length; i++) {
            slib.assert(n == arguments[1].length);
        }
        for (var i = 0; i < n; i++) {
            var items = [];
            for (var j = 0; j < n; j++) {
                items.push(arguments[j][i]);
            }
            result.push(items);
        }
    }
    return result;
}


// preload images
slib.preloadImages = function(urls)
{
    var count = 0;
    urls.forEach(function(url) {
        $('<img />').one('load', function() {
            count = count + 1;
        })
        .attr('src', url)
        .each(function() {
            if (this.complete) 
                $(this).trigger('load');
        });
    });
}


// search in an array
slib.findFirst = function(arr, func)
{
     for (var i = 0; i < arr.length; i++) {
         if (func.call(null, arr[i])) {
             break;
         }
     }
     return i;
}

slib.findLast = function(arr, func)
{
    for (var i = arr.length; i >= 0; i--) {
        if (func.call(null, arr[i])) {
            break;
        }
    }
    return i;
}

slib.min = function(arr)
{
    return Math.min.apply(null, arr);
}

slib.max = function(arr)
{
    return Math.max.apply(null, arr);
}

slib.minIndex = function(arr)
{
    return arr.reduce(function(prev, cur, idx) {
        return cur < arr[prev] ? idx : prev;
    }, 0);
}

slib.maxIndex = function(arr)
{
    return arr.reduce(function(prev, cur, idx) {
        return cur > arr[prev] ? idx : prev;
    }, 0);
}


// sorting
slib.invertIndices = function(arr)
{
    result = [];
    result.length = arr.length;
    for (var i = 0; i < result.length; i++) {
        result[arr[i]] = i;
    }
    return result;
}

slib.sortIndices = function(arr)
{
    function compare(a, b)
    {
        if (Array.isArray(a)) {
            for (var i = 0; i < a.length; i++) {
                var result = compare(a[i], b[i]);
                if (result != 0) {
                    return result;
                }
            }
            return 0;
        } else {
            if (a == b) {
                return 0;
            } else {
                return a > b ? 1 : -1;
            }
        }
    }

    return arr.map(function(x, i) { return [x, i]; })
            .sort(function(a, b) { return compare(a, b); })
            .map(function(x) { return x[1]; });
}

// element fetching
slib.fetchElems = function(arr, indices)
{
    var result = [];
    result.length = indices.length;
    for (var i = 0; i < indices.length; i++) {
        result[i] = arr[indices[i]];
    }
    return result;
}

// function argument binding
slib.bind = function(indices, args)
{
    // TODO: complete this function
}

// misc functions
slib.gotoURL = function(url)
{
    window.location.href = url;
}

slib.removeEmptyElems = function(arr)
{
    return arr.filter(function(x) { return x.length > 0; });
}

