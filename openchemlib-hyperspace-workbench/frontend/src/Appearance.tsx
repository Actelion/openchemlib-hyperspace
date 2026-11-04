import {createContext,useContext,useLayoutEffect,useState} from 'react';
import type {ReactNode} from 'react';

export type Appearance='light'|'dark';
const AppearanceContext=createContext<{theme:Appearance;toggleTheme:()=>void}>({theme:'light',toggleTheme:()=>{}});
function initialAppearance():Appearance {
  try {const stored=localStorage.getItem('hyperspace-theme');if(stored==='light'||stored==='dark')return stored;}catch{/* Storage can be unavailable in a private browser. */}
  return window.matchMedia('(prefers-color-scheme: dark)').matches?'dark':'light';
}
export function AppearanceProvider({children}:{children:ReactNode}) {
  const [theme,setTheme]=useState<Appearance>(initialAppearance);
  useLayoutEffect(()=>{
    document.documentElement.dataset.theme=theme;
    document.querySelector('meta[name="theme-color"]')?.setAttribute('content',theme==='dark'?'#111f28':'#ffffff');
    try {localStorage.setItem('hyperspace-theme',theme);}catch{/* Appearance still works without persistence. */}
  },[theme]);
  return <AppearanceContext.Provider value={{theme,toggleTheme:()=>setTheme(current=>current==='light'?'dark':'light')}}>{children}</AppearanceContext.Provider>;
}
export const useAppearance=()=>useContext(AppearanceContext);
