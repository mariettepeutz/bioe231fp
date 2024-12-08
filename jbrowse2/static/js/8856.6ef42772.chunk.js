"use strict";(globalThis.webpackChunk_jbrowse_web=globalThis.webpackChunk_jbrowse_web||[]).push([[8856],{28856:(e,t,i)=>{i.d(t,{doAfterAttach:()=>f});var n=i(42489),o=i(36422),a=i(99546),s=i(95095),r=i(42781),l=i(30385);function f(e){(0,o.addDisposer)(e,(0,n.autorun)((()=>{const t=(0,a.getContainingView)(e);if(!t.initialized)return;const i=e.mainCanvas?.getContext("2d"),n=e.cigarClickMapCanvas?.getContext("2d");if(!i||!n)return;const o=e.height,s=t.width;i.clearRect(0,0,s,o),n.clearRect(0,0,s,o),(0,l.Ww)(e,i,n)}))),(0,o.addDisposer)(e,(0,n.autorun)((()=>{(0,a.getContainingView)(e).initialized&&(0,l.C4)(e)}))),(0,o.addDisposer)(e,(0,n.reaction)((()=>{const t=(0,a.getContainingView)(e);return{bpPerPx:t.views.map((e=>e.bpPerPx)),displayedRegions:JSON.stringify(t.views.map((e=>e.displayedRegions))),features:e.features,initialized:t.initialized}}),(({initialized:t})=>{if(!t)return;const{level:i}=e,{assemblyManager:n}=(0,a.getSession)(e),l=(0,a.getContainingView)(e).views.map((e=>({...(0,o.getSnapshot)(e),width:e.width,staticBlocks:e.staticBlocks,interRegionPaddingWidth:e.interRegionPaddingWidth,minimumBlockWidth:e.minimumBlockWidth}))),f=[],c=e.features||[];for(const e of c){const t=e.get("mate");let o=e.get("start"),a=e.get("end");const c=t.start,g=t.end;-1===e.get("strand")&&([a,o]=[o,a]);const d=n.get(e.get("assemblyName")),h=n.get(t.assemblyName),u=e.get("refName"),m=t.refName,v=d?.getCanonicalRefName(u)||u,C=h?.getCanonicalRefName(m)||m,b=l[i],p=l[i+1],w=(0,s.eB)({self:b,refName:v,coord:o}),M=(0,s.eB)({self:b,refName:v,coord:a}),x=(0,s.eB)({self:p,refName:C,coord:c}),P=(0,s.eB)({self:p,refName:C,coord:g});if(void 0===w||void 0===M||void 0===x||void 0===P)continue;const k=e.get("CIGAR");f.push({p11:w,p12:M,p21:x,p22:P,f:e,cigar:r.aF.parseCigar(k)})}e.setFeatPositions(f)}),{fireImmediately:!0}))}},79610:(e,t,i)=>{i.d(t,{$2:()=>s,Eg:()=>f,WT:()=>r,f0:()=>l,mr:()=>a});var n=i(99546),o=i(30385);function a({feature:e,ctx:t,offsets:i,level:o,cb:a,height:r,drawCurves:l,oobLimit:f,viewWidth:c,hideTiny:g}){const{p11:d,p12:h,p21:u,p22:m}=e,v=d.offsetPx-i[o],C=h.offsetPx-i[o],b=u.offsetPx-i[o+1],p=m.offsetPx-i[o+1],w=Math.abs(C-v),M=Math.abs(p-b),x=r,P=(x-0)/2,k=Math.min(b,p),y=Math.max(b,p);(0,n.doesIntersect2)(k,y,-f,c+f)&&(w<=1&&M<=1?g||(t.beginPath(),t.moveTo(v,0),l?t.bezierCurveTo(v,P,b,P,b,x):t.lineTo(b,x),t.stroke()):(s(t,v,C,0,p,b,x,P,l),a(t)))}function s(e,t,i,n,o,a,s,r,l){l?function(e,t,i,n,o,a,s,r){const l=Math.abs(t-i),f=Math.abs(t-i);if(l<5&&f<5&&i<t&&Math.abs(t-o)>100){const e=t;t=i,i=e}e.beginPath(),e.moveTo(t,n),e.lineTo(i,n),e.bezierCurveTo(i,r,o,r,o,s),e.lineTo(a,s),e.bezierCurveTo(a,r,t,r,t,n),e.closePath()}(e,t,i,n,o,a,s,r):function(e,t,i,n,o,a,s){e.beginPath(),e.moveTo(t,n),e.lineTo(i,n),e.lineTo(o,s),e.lineTo(a,s),e.closePath()}(e,t,i,n,o,a,s)}function r(e,t){const i=(0,n.getContainingView)(t),a=(0,n.getContainingTrack)(t),{featPositions:s,numFeats:r,clickMapCanvas:l,cigarClickMapCanvas:f,level:c}=t;if(!l||!f)return;const g=l.getBoundingClientRect(),d=l.getContext("2d"),h=f.getContext("2d");if(!d||!h)return;const u=e.clientX-g.left,m=e.clientY-g.top,[v,C,b]=d.getImageData(u,m,1,1).data,p=Math.floor(o.xx/r),w=s[(0,o.OX)(v,C,b,p)];if(w){const{f:e}=w;t.setClickId(e.id());const o=(0,n.getSession)(t);(0,n.isSessionModelWithWidgets)(o)&&o.showWidget(o.addWidget("SyntenyFeatureWidget","syntenyFeature",{view:i,track:a,featureData:e.toJSON(),level:c}))}return w}function l(e,t,i){e.preventDefault();const n=t.clickMapCanvas,a=t.cigarClickMapCanvas;if(!n||!a)return;const s=n.getBoundingClientRect(),r=n.getContext("2d"),l=a.getContext("2d");if(!r||!l)return;const{clientX:f,clientY:c}=e,g=f-s.left,d=c-s.top,[h,u,m]=r.getImageData(g,d,1,1).data,v=Math.floor(o.xx/t.numFeats),C=(0,o.OX)(h,u,m,v),b=t.featPositions[C];b&&(t.setClickId(b.f.id()),i({clientX:f,clientY:c,feature:b}))}function f({feature:e,cigarOp:t,cigarOpLen:i}){const o=e.toJSON(),a=o.mate,s=o.end-o.start,r=a.end-a.start,l=o.identity,f=o.name,c=a.name;return[`Loc1: ${(0,n.assembleLocString)(o)}`,`Loc2: ${(0,n.assembleLocString)(a)}`,`Inverted: ${-1===o.strand}`,`Query len: ${s.toLocaleString("en-US")}`,`Target len: ${r.toLocaleString("en-US")}`,l?`Identity: ${l.toPrecision(2)}`:"",t?`CIGAR operator: ${t}${i}`:"",f?`Name 1: ${f}`:"",c?`Name 1: ${c}`:""].filter((e=>!!e)).join("<br/>")}},30385:(e,t,i)=>{i.d(t,{C4:()=>d,OX:()=>c,Ww:()=>g,xx:()=>a});var n=i(99546),o=i(79610);const a=16581375;function s(e){return`rgb(${Math.floor(e/65025)%255},${Math.floor(e/255)%255},${e%255})`}const r={I:"#ff03",N:"#0a03",D:"#00f3",X:"brown",M:"#f003","=":"#f003"},l=3,f=1600;function c(e,t,i,n){return Math.floor((255*e*255+255*t+i-1)/n)}function g(e,t,i){const c=(0,n.getContainingView)(e),g=c.drawCurves,d=c.drawCIGAR,{level:h,height:u,featPositions:m}=e,v=c.width,C=c.views.map((e=>e.bpPerPx));i&&(i.imageSmoothingEnabled=!1),t.beginPath();const b=c.views.map((e=>e.offsetPx)),p=Math.floor(a/m.length);t.fillStyle=r.M,t.strokeStyle=r.M;for(const{p11:e,p12:i,p21:n,p22:o}of m){const a=e.offsetPx-b[h],s=i.offsetPx-b[h],r=n.offsetPx-b[h+1],c=o.offsetPx-b[h+1],d=Math.abs(s-a),m=Math.abs(c-r),C=0,p=u,w=(p-C)/2;d<=l&&m<=l&&r<v+f&&r>-f&&(t.moveTo(a,C),g?t.bezierCurveTo(a,w,r,w,r,p):t.lineTo(r,p))}t.stroke(),t.fillStyle=r.M,t.strokeStyle=r.M;for(const{p11:e,p12:p,p21:w,p22:M,f:x,cigar:P}of m){const m=e.offsetPx-b[h],k=p.offsetPx-b[h],y=w.offsetPx-b[h+1],S=M.offsetPx-b[h+1],T=Math.abs(k-m),N=Math.abs(S-y),$=Math.min(y,S),I=Math.max(y,S),W=0,R=u,L=(R-W)/2;if(!(T<=l&&N<=l)&&(0,n.doesIntersect2)($,I,-f,c.width+f)){const e=x.get("strand"),n=-1===e?k:m,l=n<(-1===e?m:k)?1:-1,f=(y<S?1:-1)*e;let c=n,u=-1===e?S:y;if(P.length&&d){let e=!1,n=0,d=0;const m=Math.floor(a/P.length);for(let a=0;a<P.length;a+=2){const b=a*m+1,p=+P[a],w=P[a+1];e||(n=c,d=u);const M=p/C[h],x=p/C[h+1];if("M"===w||"="===w||"X"===w?(c+=M*l,u+=x*f):"D"===w||"N"===w?c+=M*l:"I"===w&&(u+=x*f),!(Math.max(n,d,c,u)<0||Math.min(n,d,c,u)>v)){const l=a<P.length-2;Math.abs(c-n)<=1&&Math.abs(u-d)<=1&&l?e=!0:(t.fillStyle=r[e&&M>1||x>1?w:"M"],e=!1,(0,o.$2)(t,n,c,W,u,d,R,L,g),t.fill(),i&&(i.fillStyle=s(b),(0,o.$2)(i,n,c,W,u,d,R,L,g),i.fill()))}}}else(0,o.$2)(t,m,k,W,S,y,R,L,g),t.fill()}}const w=e.clickMapCanvas?.getContext("2d");if(w){w.imageSmoothingEnabled=!1,w.clearRect(0,0,v,u);for(let e=0;e<m.length;e++){const t=m[e],i=e*p+1;w.fillStyle=s(i),(0,o.mr)({cb:e=>{e.fill()},feature:t,ctx:w,drawCurves:g,level:h,offsets:b,oobLimit:f,viewWidth:c.width,hideTiny:!0,height:u})}}}function d(e){const{level:t,clickId:i,mouseoverId:a}=e,s=(0,n.getContainingView)(e),r=s.drawCurves,l=e.height,c=s.width,g=e.mouseoverCanvas?.getContext("2d"),d=s.views.map((e=>e.offsetPx));if(!g)return;g.resetTransform(),g.scale(1,1),g.clearRect(0,0,c,l),g.strokeStyle="rgba(0, 0, 0, 0.9)",g.fillStyle="rgba(0, 0, 0, 0.1)";const h=e.featMap[a||""];h&&(0,o.mr)({cb:e=>{e.fill()},feature:h,level:t,ctx:g,oobLimit:f,viewWidth:s.width,drawCurves:r,offsets:d,height:l});const u=e.featMap[i||""];u&&(0,o.mr)({cb:e=>{e.stroke()},feature:u,ctx:g,level:t,oobLimit:f,viewWidth:s.width,drawCurves:r,offsets:d,height:l})}}}]);
//# sourceMappingURL=8856.6ef42772.chunk.js.map